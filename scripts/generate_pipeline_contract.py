#!/usr/bin/env python3
"""Pipeline Contract & API Surface Generator.

Performs static AST analysis on src/ to generate docs/PIPELINE_AI_CONTRACT.md
and docs/pipeline_api_surface.json. Serves as the immutable reference for all
pipeline engines, exact signatures, and orchestration invariants.
"""

import ast
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
DOCS_DIR = PROJECT_ROOT / "docs"


def get_node_source(node: Optional[ast.AST]) -> str:
    """Safely converts an AST type annotation node to its string representation."""
    if node is None:
        return "Any"
    if hasattr(ast, "unparse"):
        return ast.unparse(node)
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return f"{get_node_source(node.value)}.{node.attr}"
    if isinstance(node, ast.Subscript):
        return f"{get_node_source(node.value)}[{get_node_source(node.slice)}]"
    if isinstance(node, ast.Constant):
        return repr(node.value)
    return "..."


def parse_function_signature(fn_node: ast.FunctionDef) -> Dict[str, Any]:
    """Extracts exact parameters, default values, type annotations, and docstring."""
    decorators = [get_node_source(d) for d in fn_node.decorator_list]
    args_info = []

    num_defaults = len(fn_node.args.defaults)
    pos_args = fn_node.args.args
    defaults_offset = len(pos_args) - num_defaults

    for idx, arg in enumerate(pos_args):
        default_val = None
        if idx >= defaults_offset:
            default_val = get_node_source(fn_node.args.defaults[idx - defaults_offset])
        args_info.append({
            "name": arg.arg,
            "type": get_node_source(arg.annotation) if arg.annotation else "Any",
            "default": default_val,
        })

    if fn_node.args.vararg:
        args_info.append({"name": f"*{fn_node.args.vararg.arg}", "type": "Any", "default": None})

    for idx, arg in enumerate(fn_node.args.kwonlyargs):
        default_val = None
        if idx < len(fn_node.args.kw_defaults) and fn_node.args.kw_defaults[idx]:
            default_val = get_node_source(fn_node.args.kw_defaults[idx])
        args_info.append({
            "name": arg.arg,
            "type": get_node_source(arg.annotation) if arg.annotation else "Any",
            "default": default_val,
        })

    if fn_node.args.kwarg:
        args_info.append({"name": f"**{fn_node.args.kwarg.arg}", "type": "Any", "default": None})

    return_type = get_node_source(fn_node.returns) if fn_node.returns else "None"
    docstring = ast.get_docstring(fn_node) or ""
    summary_doc = docstring.strip().split("\n")[0] if docstring else ""

    return {
        "name": fn_node.name,
        "decorators": decorators,
        "is_classmethod": "classmethod" in decorators,
        "is_staticmethod": "staticmethod" in decorators,
        "parameters": args_info,
        "return_type": return_type,
        "docstring_summary": summary_doc,
    }


def parse_class_definition(cls_node: ast.ClassDef) -> Dict[str, Any]:
    """Extracts class metadata, bases, and public/private methods."""
    bases = [get_node_source(b) for b in cls_node.bases]
    docstring = ast.get_docstring(cls_node) or ""
    summary_doc = docstring.strip().split("\n")[0] if docstring else ""

    methods = []
    for item in cls_node.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            methods.append(parse_function_signature(item))

    return {
        "name": cls_node.name,
        "bases": bases,
        "docstring_summary": summary_doc,
        "methods": methods,
    }


def extract_module_ast(file_path: Path) -> Dict[str, Any]:
    """Extracts all classes and top-level functions from a Python source file."""
    rel_path = file_path.relative_to(PROJECT_ROOT)
    try:
        content = file_path.read_text(encoding="utf-8")
        tree = ast.parse(content, filename=str(file_path))
    except Exception as exc:
        return {"file": str(rel_path), "error": str(exc), "classes": [], "functions": []}

    classes = []
    functions = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            classes.append(parse_class_definition(node))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.append(parse_function_signature(node))

    return {
        "file": str(rel_path),
        "module_path": ".".join(rel_path.with_suffix("").parts),
        "classes": classes,
        "functions": functions,
    }


def generate_contract_markdown(modules_data: List[Dict[str, Any]]) -> str:
    """Builds the comprehensive markdown contract file."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    lines = [
        "# COMPLETE PIPELINE ARCHITECTURE, TEST CONTRACT & API SURFACE",
        f"> **Generated:** {timestamp} | **Root:** `amtdna-variant-pipeline`",
        "",
        "## 1. BLOCK-BY-BLOCK CONTRACT & INVARIANTS",
        "",
        "### Block 1: Reference Genome Management",
        "- **Engine**: `src.modules.alignment.reference_manager.ReferenceManager`",
        "- **Primary Entrypoint**: `prepare_and_validate_reference(fasta_path, expected_contig, strict_contig, expected_length, allow_n, auto_index, log_dir)`",
        "- **Invariants**:",
        "  1. `allow_n=True` is strictly required for rCRS (NC_012920.1) due to position 3107.",
        "  2. Validates physical FAI and 5 BWA indices (`.bwt`, `.pac`, `.ann`, `.amb`, `.sa`).",
        "  3. Emits validated `ReferenceValidationResult` containing canonical `ReferenceInfo`.",
        "",
        "### Block 2: FASTQ Discovery, Integrity & Preprocessing",
        "- **Engines**:",
        "  - Discovery: `src.modules.preprocessing.discovery.SampleDiscoveryEngine`",
        "  - Integrity: `src.modules.preprocessing.integrity.IntegrityEngine`",
        "  - Preprocessing: `src.modules.preprocessing.trimmer.PreprocessingEngine`",
        "- **Primary Entrypoints**:",
        "  - `SampleDiscoveryEngine.discover_files(input_dir)`",
        "  - `SampleDiscoveryEngine.pair_and_validate_fastq(discovered_files, allow_missing_mates, single_end_mode)`",
        "  - `SampleDiscoveryEngine.build_sample_records(paired_units, orphan_files)`",
        "  - `IntegrityEngine.verify_sample_integrity(sample)`",
        "  - `PreprocessingEngine.process_sample(sample, config)`",
        "- **Invariants**:",
        "  1. AdapterRemoval trims adapters, merges overlapping reads, and isolates branches (`collapsed`, `pair1`/`pair2`, `single`).",
        "  2. Verifies SHA-256 fingerprint on raw inputs.",
        "",
        "### Block 3: Alignment, Deduplication, Filtering & QC",
        "- **Orchestrator**: `src.modules.alignment.orchestrator.AlignmentOrchestrator`",
        "- **Primary Entrypoint**: `run_sample_alignment(sample, reference_info, output_dir, log_dir, policy, threads, ...)`",
        "- **Invariants**:",
        "  1. Multi-branch alignment with distinct Read Groups (`@RG LB:PAIRED`, `@RG LB:MERGED`, `@RG LB:SINGLE`).",
        "  2. Canonical ordering: Merge branches -> Samtools collate/fixmate -m -> Coordinate sort -> Markdup -> MAPQ filtering.",
        "  3. Outputs `filtered.dedup.bam` (MAPQ >= 30, flags 0xF04 filtered) and indexes.",
        "",
        "### Block 4: aDNA Damage Authentication & Profiling",
        "- **Orchestrator**: `src.modules.damage.orchestrator.DamageOrchestrator`",
        "- **Primary Entrypoint**: `run_sample_damage_pipeline(sample_id, reference_info, output_dir, log_dir, policy, dedup_bam, ...)`",
        "- **Invariants**:",
        "  1. Executes mapDamage2 statistical quantification (5' C->T and 3' G->A).",
        "  2. Generates base-quality rescaled BAM (`rescaled.bam`) and PMD score metrics.",
        "  3. Downstream Block 5 MUST receive the unrescaled deduplicated BAM.",
        "",
        "### Block 5: Two-Tier Contamination Triage & Bayesian Modeling",
        "- **Batch Runner**: `src.modules.contamination.batch_runner.ContaminationBatchRunner`",
        "- **Primary Entrypoint**: `run_cohort_pipeline(discovered_samples, reference_fasta, reference_contig, reference_length_bp, output_dir, config, resource_alloc)`",
        "- **Invariants**:",
        "  1. FastContaminationScreen evaluates 354 diagnostic sites concurrently.",
        "  2. `FAST_CLEAR` skips Schmutzi Bayesian modeling completely.",
        "  3. `ESCALATE` and `INSUFFICIENT_DATA` queue bounded execution of Schmutzi in isolated containers.",
        "  4. Global host CPU budget fraction is strictly enforced across all steps.",
        "",
        "---",
        "",
        "## 2. STRICT OPERATIONAL RULES (DO NOT LIST)",
        "1. **DO NOT** instantiate ad-hoc command wrappers; invoke official domain classes directly.",
        "2. **DO NOT** pass `reference_info` to `prepare_and_validate_reference(fasta_path=...)`.",
        "3. **DO NOT** pass rescaled BAMs to Schmutzi in Block 5.",
        "4. **DO NOT** bypass Schmutzi if FastScreen returns `INSUFFICIENT_DATA`.",
        "5. **DO NOT** open `.sai` alignment output as text streams (handle raw binary only).",
        "6. **DO NOT** mark a stage as SUCCESS without verifying physical file existence and non-zero size on disk.",
        "",
        "---",
        "",
        "## 3. EXACT CODEBASE API SURFACE (CLASSES, METHODS & SIGNATURES)",
        "",
    ]

    for mod in modules_data:
        if mod.get("error"):
            continue
        if not mod["classes"] and not mod["functions"]:
            continue

        lines.append(f"### Module `{mod['module_path']}`")
        lines.append(f"**File**: `{mod['file']}`\n")

        for cls in mod["classes"]:
            bases_str = f"({', '.join(cls['bases'])})" if cls["bases"] else ""
            lines.append(f"#### `class {cls['name']}{bases_str}`")
            if cls["docstring_summary"]:
                lines.append(f"> {cls['docstring_summary']}\n")

            for m in cls["methods"]:
                # Filter out private methods except __init__ and __post_init__
                if m["name"].startswith("_") and m["name"] not in ("__init__", "__post_init__"):
                    continue

                params_str = []
                for p in m["parameters"]:
                    p_str = p["name"]
                    if p["type"] != "Any":
                        p_str += f": {p['type']}"
                    if p["default"] is not None:
                        p_str += f" = {p['default']}"
                    params_str.append(p_str)

                dec_prefix = ""
                if m["is_classmethod"]:
                    dec_prefix = "@classmethod\n"
                elif m["is_staticmethod"]:
                    dec_prefix = "@staticmethod\n"

                lines.append("```python")
                lines.append(f"{dec_prefix}def {m['name']}({', '.join(params_str)}) -> {m['return_type']}:")
                if m["docstring_summary"]:
                    lines.append(f'    """{m["docstring_summary"]}"""')
                lines.append("```\n")

        for fn in mod["functions"]:
            if fn["name"].startswith("_"):
                continue
            params_str = []
            for p in fn["parameters"]:
                p_str = p["name"]
                if p["type"] != "Any":
                    p_str += f": {p['type']}"
                if p["default"] is not None:
                    p_str += f" = {p['default']}"
                params_str.append(p_str)

            lines.append("```python")
            lines.append(f"def {fn['name']}({', '.join(params_str)}) -> {fn['return_type']}:")
            if fn["docstring_summary"]:
                lines.append(f'    """{fn["docstring_summary"]}"""')
            lines.append("```\n")

        lines.append("---\n")

    return "\n".join(lines)


def main() -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    py_files = sorted(list(SRC_DIR.rglob("*.py")))
    py_files = [f for f in py_files if "__pycache__" not in f.parts and not f.name.endswith(".save")]

    print(f"[*] Extracting AST API surface from {len(py_files)} Python source files...")
    all_modules = []
    for f in py_files:
        mod_data = extract_module_ast(f)
        all_modules.append(mod_data)

    # 1. Export JSON surface
    json_path = DOCS_DIR / "pipeline_api_surface.json"
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump({"generated_at": datetime.now(timezone.utc).isoformat(), "modules": all_modules}, jf, indent=2)
    print(f"✅ Generated API JSON: {json_path}")

    # 2. Export Markdown Contract
    md_path = DOCS_DIR / "PIPELINE_AI_CONTRACT.md"
    md_content = generate_contract_markdown(all_modules)
    md_path.write_text(md_content, encoding="utf-8")
    print(f"✅ Generated AI Contract: {md_path}")
    print(f"[*] Total lines generated: {len(md_content.splitlines())}")


if __name__ == "__main__":
    main()