"""
Forensic mtDNA Pipeline - Production Graphical HUD.
Dispatches background jobs, manages atomic lifecycle states,
synchronizes BAM/FASTQ path defaults dynamically via clean callbacks,
provides interactive pre-flight cache resolution (Reuse / Overwrite / Cancel),
provides adaptive resilient TSV reads with memory fallback,
and loads multi-palette themes (gradients/images) from config/themes.yaml.
"""

from collections.abc import Hashable
import json
import os
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import streamlit as st
import yaml

CONFIG_FILE = PROJECT_ROOT / "config" / "default_config.yaml"
THEMES_FILE = PROJECT_ROOT / "config" / "themes.yaml"

from src.gui.log_highlighter import LogHighlighter
from src.gui.process_manager import (
    RUN_LOG_FILE,
    PathValidator,
    ProcessManager,
    StopJobStatus,
)
from src.pipeline.stage_tracker import StageTracker

# ------------------------------------------------------------------------------
# Dynamic Theme Loader from config/themes.yaml
# ------------------------------------------------------------------------------
DEFAULT_THEMES = {
    "Midnight Petrol (Bio-Abyss Gradient)": {
        "bg_style": "linear-gradient(135deg, #09121D 0%, #0D1B2A 50%, #162638 100%)",
        "sidebar_bg": "#08121E",
        "card_bg": "rgba(19, 39, 61, 0.75)",
        "text_color": "#E0F2FE",
        "muted_text": "#7DD3FC",
        "code_bg": "#04090F",
        "code_text": "#FFFFFF",
        "accent": "#06B6D4",
    }
}

def load_themes() -> Dict[str, Dict[str, str]]:
    if THEMES_FILE.is_file():
        try:
            with open(THEMES_FILE, "r", encoding="utf-8") as fp:
                data = yaml.safe_load(fp) or {}
                loaded = data.get("themes", {})
                if isinstance(loaded, dict) and loaded:
                    return loaded
        except Exception:
            pass
    return DEFAULT_THEMES

THEMES = load_themes()

# ------------------------------------------------------------------------------
# Resilient TSV Reader with Adaptive Backoff & Memory Fallback
# ------------------------------------------------------------------------------
def read_tsv_resilient(
    filepath: Path,
    retries: int = 8,
    initial_delay: float = 0.05,
    max_delay: float = 0.40,
    **kwargs: Any,
) -> Optional[pd.DataFrame]:
    cache: dict[str, pd.DataFrame] = st.session_state.setdefault("_tsv_memory_cache", {})
    cache_key = str(filepath.resolve()) if filepath.is_file() else str(filepath)

    delay = initial_delay
    for _ in range(retries):
        try:
            if not filepath.is_file():
                time.sleep(delay)
                delay = min(delay * 1.5, max_delay)
                continue
            df = pd.read_csv(filepath, sep="\t", **kwargs)
            cache[cache_key] = df
            return df
        except (FileNotFoundError, PermissionError, pd.errors.EmptyDataError, pd.errors.ParserError):
            time.sleep(delay)
            delay = min(delay * 1.5, max_delay)

    return cache.get(cache_key)


# Load output directories and default cache mode from configuration
config_data = {}
if CONFIG_FILE.is_file():
    with open(CONFIG_FILE, "r", encoding="utf-8") as fp:
        config_data = yaml.safe_load(fp) or {}

out_cfg = config_data.get("output", {})
results_dir = PROJECT_ROOT / out_cfg.get("results_dir", "results")
operational_files_dir = PROJECT_ROOT / out_cfg.get("operational_files_dir", "operational_files")
default_yaml_cache_mode = str(out_cfg.get("cache_mode", "prompt")).lower()

master_dir = results_dir / "master_summary"
reports_dir = results_dir / "reports_markdown"
mutations_dir = results_dir / "mutations"

# ------------------------------------------------------------------------------
# Interactive Modal Dialog for Pre-Flight Cache Decision
# ------------------------------------------------------------------------------
@st.dialog("⚠️ Existing Cache Detected")
def prompt_cache_decision_dialog(
    detected_cache: Dict[str, List[str]],
    data_type: str,
    input_dir_arg: str,
    cpu_fraction: float,
    log_level: str,
) -> None:
    st.warning(
        f"Found pre-existing cached stages for **{len(detected_cache)} sample(s)** "
        f"in `{operational_files_dir.name}/`."
    )

    with st.expander("🔍 View Cached Samples & Stages", expanded=True):
        for s_id, stages in detected_cache.items():
            st.markdown(f"- **`{s_id}`**: `{', '.join(stages)}`")

    st.markdown("**How would you like to proceed?**")

    col_reuse, col_overwrite, col_cancel = st.columns(3)

    with col_reuse:
        if st.button("♻️ Reuse Cache", type="primary", width="stretch"):
            success, msg = ProcessManager.start_job(
                data_type=data_type,
                input_dir_arg=input_dir_arg,
                cpu_fraction=cpu_fraction,
                log_level=log_level,
                cache_mode="reuse",
            )
            if success:
                st.session_state["_launch_feedback"] = ("success", msg)
            else:
                st.session_state["_launch_feedback"] = ("error", msg)
            st.rerun()

    with col_overwrite:
        if st.button("🔥 Overwrite Cache", type="secondary", width="stretch"):
            success, msg = ProcessManager.start_job(
                data_type=data_type,
                input_dir_arg=input_dir_arg,
                cpu_fraction=cpu_fraction,
                log_level=log_level,
                cache_mode="overwrite",
            )
            if success:
                st.session_state["_launch_feedback"] = ("success", msg)
            else:
                st.session_state["_launch_feedback"] = ("error", msg)
            st.rerun()

    with col_cancel:
        if st.button("🛑 Stop / Cancel", width="stretch"):
            st.session_state["_launch_feedback"] = ("info", "Job launch cancelled by user.")
            st.rerun()


# ------------------------------------------------------------------------------
# Sidebar Controls & Theme Selector
# ------------------------------------------------------------------------------
is_running, active_job = ProcessManager.is_running()

st.sidebar.title("HUD Configuration")

theme_choice = st.sidebar.selectbox(
    "🎨 Theme Palette",
    list(THEMES.keys()),
    index=0,
    help="Loaded dynamically from config/themes.yaml",
)
t = THEMES[theme_choice]

# Support both flat colors and full gradients/images
bg_css_value = t.get("bg_style") or t.get("bg_color", "#0B1120")

# Inject Scoped CSS with Glassmorphism for cards
st.markdown(
    f"""
    <style>
    /* Main Background (Gradients & Images supported) */
    .stApp {{
        background: {bg_css_value} !important;
        background-attachment: fixed !important;
        background-size: cover !important;
        color: {t['text_color']} !important;
    }}
    .stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp p, .stApp span {{
        color: {t['text_color']};
    }}
    
    /* Sidebar */
    [data-testid="stSidebar"] {{
        background-color: {t['sidebar_bg']} !important;
        border-right: 1px solid {t['accent']}33 !important;
    }}
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3,
    [data-testid="stSidebar"] label, [data-testid="stSidebar"] p {{
        color: {t['text_color']} !important;
    }}

    /* Tabs Styling */
    button[data-baseweb="tab"] {{
        background-color: transparent !important;
        color: {t['muted_text']} !important;
        font-weight: 600 !important;
    }}
    button[data-baseweb="tab"][aria-selected="true"] {{
        color: {t['accent']} !important;
        border-bottom-color: {t['accent']} !important;
    }}

    /* Metric Cards with Backdrop Glass Effect */
    div[data-testid="stMetric"] {{
        background: {t['card_bg']} !important;
        backdrop-filter: blur(10px) !important;
        -webkit-backdrop-filter: blur(10px) !important;
        border: 1px solid {t['accent']}44 !important;
        padding: 12px 18px !important;
        border-radius: 8px !important;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15) !important;
    }}
    div[data-testid="stMetricValue"] {{
        color: {t['accent']} !important;
    }}
    </style>
    """,
    unsafe_allow_html=True,
)

st.sidebar.markdown("---")
st.sidebar.subheader("Pipeline Controls")

if is_running:
    st.sidebar.warning(f"🔒 Job is running (PID: {active_job.pid if active_job else 'active'}). Controls locked.")

default_bam_dir = str(PROJECT_ROOT / "data" / "raw_bam")
default_fastq_dir = str(PROJECT_ROOT / "data" / "raw_fastq")

if "input_dir_text" not in st.session_state:
    st.session_state["input_dir_text"] = default_bam_dir

def on_data_type_change() -> None:
    new_mode = st.session_state["data_type_selector"]
    st.session_state["input_dir_text"] = default_bam_dir if new_mode == "BAM" else default_fastq_dir

data_type_choice = st.sidebar.radio(
    "Sequence Input Mode",
    ["BAM", "FASTQ"],
    index=0,
    disabled=is_running,
    key="data_type_selector",
    on_change=on_data_type_change,
)

raw_input_dir = st.sidebar.text_input(
    "Input Directory (Windows or WSL)",
    key="input_dir_text",
    disabled=is_running,
)

dir_exists, runtime_target_dir, pipeline_posix_arg = PathValidator.resolve_paths(raw_input_dir)

cpu_percent = st.sidebar.slider(
    "Host CPU Budget",
    min_value=10,
    max_value=100,
    value=50,
    step=5,
    format="%d%%",
    disabled=is_running,
)
cpu_fraction = cpu_percent / 100.0

cache_mode_options = ["prompt", "reuse", "overwrite"]
default_cache_idx = (
    cache_mode_options.index(default_yaml_cache_mode)
    if default_yaml_cache_mode in cache_mode_options
    else 0
)
cache_mode_labels = {
    "prompt": "Ask on Start (Prompt)",
    "reuse": "Reuse Existing Cache",
    "overwrite": "Overwrite (Purge & Recompute)",
}

cache_policy_choice = st.sidebar.selectbox(
    "Cache Handling Policy",
    options=cache_mode_options,
    index=default_cache_idx,
    format_func=lambda k: cache_mode_labels.get(k, k),
    disabled=is_running,
    help="Default loaded from config/default_config.yaml (output.cache_mode).",
)

log_level = st.sidebar.selectbox(
    "Terminal Log Verbosity",
    ["INFO", "DEBUG", "WARNING", "ERROR", "CRITICAL"],
    index=0,
    disabled=is_running,
)

st.sidebar.markdown("---")
st.sidebar.caption(f"Pipeline Arg: `{pipeline_posix_arg}`")

# ------------------------------------------------------------------------------
# Main Dashboard
# ------------------------------------------------------------------------------
st.title("🧬 mtDNA Forensic & Ancient DNA Pipeline HUD")

tab_run, tab_cohort, tab_sample, tab_glossary = st.tabs([
    "🚀 Execution & Live Monitor",
    "📊 Cohort QC Dashboard",
    "🔍 Sample Deep-Dive",
    "📖 Diagnostic Glossary",
])

# ------------------------------------------------------------------------------
# Fragment: Isolated Real-time Log Streamer
# ------------------------------------------------------------------------------
@st.fragment(run_every="3s")
def render_live_log(was_active_at_mount: bool, active_theme: dict) -> None:
    running, _ = ProcessManager.is_running()

    if was_active_at_mount and not running:
        last_c = ProcessManager.get_last_exit_code()
        if last_c == 0:
            st.toast("🎉 Pipeline execution finished successfully!", icon="✅")
        else:
            st.toast(f"⚠️ Pipeline finished with exit code {last_c}", icon="🔴")

    if running:
        st.info("⚡ Pipeline is actively running in background. You can safely browse tabs.")
    else:
        last_code = ProcessManager.get_last_exit_code()
        if last_code is not None:
            if last_code == 0:
                st.success("Last execution completed successfully (Exit Code: 0).")
            else:
                st.error(f"Last execution ended with errors (Exit Code: {last_code}).")

    st.markdown("### Execution Log Stream")
    lines = LogHighlighter.tail_lines_binary(RUN_LOG_FILE, n_lines=60)
    if lines:
        st.html(LogHighlighter.render_terminal_html(lines, theme=active_theme))
    else:
        st.caption("No log records available.")


# ------------------------------------------------------------------------------
# TAB 1: Execution & Process Controls
# ------------------------------------------------------------------------------
with tab_run:
    st.subheader("Process Management")

    # Display feedback from modal dialog actions if present
    if "_launch_feedback" in st.session_state:
        fb_type, fb_msg = st.session_state.pop("_launch_feedback")
        if fb_type == "success":
            st.success(fb_msg)
        elif fb_type == "error":
            st.error(fb_msg)
        else:
            st.info(fb_msg)

    detected_cache = StageTracker.inspect_operational_cache(operational_files_dir)

    c1, c2, c3 = st.columns([2, 1, 1])

    with c1:
        if dir_exists:
            pattern = "*.bam" if data_type_choice == "BAM" else "*.f*q*"
            files_found = list(runtime_target_dir.glob(pattern)) + list(runtime_target_dir.glob(pattern.upper()))
            cache_badge = (
                f" | 🗄️ {len(detected_cache)} cached sample(s) in `{operational_files_dir.name}/`"
                if detected_cache
                else " | 🆕 Clean workspace (no cache)"
            )
            st.success(f"Directory ready: {len(files_found)} candidate file(s) found{cache_badge}.")
        else:
            st.error(f"Target directory does not exist on host: `{runtime_target_dir}`")

    with c2:
        if not is_running:
            if st.button("▶ Start Pipeline", type="primary", disabled=not dir_exists):
                if cache_policy_choice == "prompt" and detected_cache:
                    prompt_cache_decision_dialog(
                        detected_cache=detected_cache,
                        data_type=data_type_choice,
                        input_dir_arg=pipeline_posix_arg,
                        cpu_fraction=cpu_fraction,
                        log_level=log_level,
                    )
                else:
                    effective_mode = "reuse" if cache_policy_choice == "prompt" else cache_policy_choice
                    success, msg = ProcessManager.start_job(
                        data_type=data_type_choice,
                        input_dir_arg=pipeline_posix_arg,
                        cpu_fraction=cpu_fraction,
                        log_level=log_level,
                        cache_mode=effective_mode,
                    )
                    if success:
                        st.success(msg)
                        st.rerun()
                    else:
                        st.error(msg)
        else:
            if st.button("🛑 Stop Pipeline", type="secondary"):
                status, msg = ProcessManager.stop_job()
                if status == StopJobStatus.STOPPED_CONFIRMED:
                    st.warning(msg)
                elif status == StopJobStatus.ALREADY_STOPPED:
                    st.info(msg)
                elif status == StopJobStatus.SECURITY_REJECTED:
                    st.error(f"🛑 {msg}")
                else:
                    st.error(msg)
                st.rerun()

    with c3:
        if st.button("🔄 Refresh Artifacts"):
            st.rerun()

    st.markdown("---")
    render_live_log(is_running, t)


# ------------------------------------------------------------------------------
# TAB 2: Cohort QC Dashboard
# ------------------------------------------------------------------------------
with tab_cohort:
    st.subheader("Cohort Quality Summary")
    summary_tsv = master_dir / "cohort_qc_summary.tsv"

    if summary_tsv.is_file():
        df_cohort = read_tsv_resilient(summary_tsv, keep_default_na=False, dtype=str)
        if df_cohort is not None and not df_cohort.empty:
            c_m1, c_m2, c_m3, c_m4 = st.columns(4)
            c_m1.metric("Total Specimens", len(df_cohort))

            status_col = df_cohort["overall_status"] if "overall_status" in df_cohort else pd.Series(dtype=str)
            c_m2.metric("Passed (PASS)", int((status_col == "PASS").sum()) if len(status_col) else "n/d")
            c_m3.metric("Warnings (WARN)", int((status_col == "WARN").sum()) if len(status_col) else "n/d")
            c_m4.metric("Failures (FAIL)", int((status_col == "FAIL").sum()) if len(status_col) else "n/d")

            st.markdown("---")
            st.dataframe(df_cohort, width="stretch")

            st.download_button(
                label="📥 Download cohort_qc_summary.tsv",
                data=summary_tsv.read_bytes(),
                file_name="cohort_qc_summary.tsv",
                mime="text/tab-separated-values",
            )
        else:
            st.warning("Cohort summary file is currently being updated by an atomic operation.")
    else:
        st.info("No cohort summary table found. Execute the pipeline to emit deliverables.")


# ------------------------------------------------------------------------------
# TAB 3: Sample Deep-Dive (Full Cohort Visibility)
# ------------------------------------------------------------------------------
with tab_sample:
    st.subheader("Individual Specimen Forensic Inspection")

    cohort_samples_map: Dict[str, Dict[str, str]] = {}
    summary_tsv = master_dir / "cohort_qc_summary.tsv"

    if summary_tsv.is_file():
        df_summary = read_tsv_resilient(summary_tsv, keep_default_na=False, dtype=str)
        if df_summary is not None and "sample_id" in df_summary.columns:
            for _, row in df_summary.iterrows():
                sid = row["sample_id"]
                cohort_samples_map[sid] = {
                    "status": row.get("overall_status", "UNKNOWN"),
                    "type": row.get("biological_type", "UNKNOWN"),
                    "haplogroup": row.get("resolved_haplogroup", row.get("primary_haplogroup", "UNKNOWN")),
                }

    if reports_dir.is_dir():
        for item in reports_dir.iterdir():
            sid = item.name if item.is_dir() else item.name.removesuffix(".report.md")
            if sid not in cohort_samples_map:
                cohort_samples_map[sid] = {
                    "status": "PASS",
                    "type": "UNKNOWN",
                    "haplogroup": "UNKNOWN",
                }

    if cohort_samples_map:
        sample_ids = sorted(list(cohort_samples_map.keys()))

        def format_sample_label(sid: Hashable) -> str:
            sid_str = str(sid)
            info = cohort_samples_map.get(sid_str, {})
            st_flag = info.get("status", "UNKNOWN")
            badge = "🟢" if st_flag == "PASS" else ("🟡" if st_flag == "WARN" else "🔴")
            return f"{badge} {sid_str} (Status: {st_flag} | HG: {info.get('haplogroup', 'UNKNOWN')})"

        selected_sample = st.selectbox(
            "Select Specimen",
            sample_ids,
            format_func=format_sample_label,
        )

        c_rep, c_mut = st.columns([3, 2])

        with c_rep:
            sample_md_file = reports_dir / selected_sample / f"{selected_sample}.report.md"
            if not sample_md_file.is_file():
                sample_md_file = reports_dir / f"{selected_sample}.report.md"

            if sample_md_file.is_file():
                md_text = sample_md_file.read_text(encoding="utf-8")
                st.markdown(md_text)
                st.download_button(
                    label=f"📥 Download {selected_sample}.report.md",
                    data=md_text.encode("utf-8"),
                    file_name=f"{selected_sample}.report.md",
                    mime="text/markdown",
                )
            else:
                sample_info = cohort_samples_map[selected_sample]
                st.error(f"### 🛑 Forensic Notice: `{selected_sample}` was Gated / Rejected")
                st.markdown(
                    f"""
                    * **Overall Status**: `{sample_info['status']}`
                    * **Biological Classification**: `{sample_info['type']}`
                    * **Reason**: This specimen was excluded from consensus generation during pipeline gating
                      (e.g., contamination exceeding critical limits or unverified sequence authenticity).
                    * *No canonical FASTA or detailed report was generated for this sample.*
                    """
                )

                single_master_json = master_dir / f"{selected_sample}.master.json"
                if single_master_json.is_file():
                    try:
                        with open(single_master_json, "r", encoding="utf-8") as jf:
                            s_data = json.load(jf)
                            st.markdown(f"**Diagnostic Narrative**: {s_data.get('narrative', 'N/A')}")
                            warns = s_data.get("warnings", [])
                            if warns:
                                st.markdown("#### Diagnostic Warnings Registered:")
                                for w in warns:
                                    st.warning(f"**[{w.get('code', 'WARNING')}]**: {w.get('message', '')}")
                    except Exception:
                        pass

        with c_mut:
            st.markdown(f"### Mutation Catalog: `{selected_sample}`")
            mut_file = mutations_dir / selected_sample / f"{selected_sample}.mutation_catalog.tsv"
            if not mut_file.is_file():
                mut_file = mutations_dir / f"{selected_sample}.mutation_catalog.tsv"

            if mut_file.is_file():
                df_mut = read_tsv_resilient(mut_file, keep_default_na=False, dtype=str)
                if df_mut is not None:
                    st.dataframe(df_mut, height=550, width="stretch")
                    st.download_button(
                        label=f"📥 Download {selected_sample}.mutation_catalog.tsv",
                        data=mut_file.read_bytes(),
                        file_name=f"{selected_sample}.mutation_catalog.tsv",
                        mime="text/tab-separated-values",
                    )
                else:
                    st.caption("Mutation catalog file is currently being updated.")
            else:
                st.caption(f"No catalog artifact produced (Sample status: {cohort_samples_map[selected_sample]['status']}).")
    else:
        st.info("No specimens detected in deliverables folder.")


# ------------------------------------------------------------------------------
# TAB 4: Diagnostic Warnings & Errors Glossary
# ------------------------------------------------------------------------------
with tab_glossary:
    st.subheader("Diagnostic Warnings & Errors Glossary")

    runtime_glossary = master_dir / "WARNINGS_GLOSSARY.md"
    static_glossary = PROJECT_ROOT / "docs" / "WARNINGS_GLOSSARY.md"

    active_glossary = None
    if runtime_glossary.is_file():
        active_glossary = runtime_glossary
    elif static_glossary.is_file():
        active_glossary = static_glossary

    if active_glossary:
        st.markdown(active_glossary.read_text(encoding="utf-8"))
    else:
        st.info(
            "Glossary documentation not found. Complete Block 8 to generate runtime documentation."
        )