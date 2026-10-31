"""
Semantic Terminal Log Tokenizer and Reverse Binary Tail Reader.
Implements bounded backwards-seeking log extraction, ANSI stripping,
and single-pass regex tokenization over tail lines.
"""

import html
import os
from pathlib import Path
import re
from typing import Dict, List, Optional, Tuple

ANSI_ESCAPE_RE = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")

# Enhanced single-pass tokenization pattern with strict exclusions
TOKEN_RE = re.compile(
    r"(?P<TIMESTAMP>\[?\b\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(?:\.\d+)?\b\]?)"
    r"|(?P<LEVEL>\[(INFO|DEBUG|WARNING|ERROR|CRITICAL)\])"
    r"|(?P<BLOCK>▶\s*BLOCK\s*\d+:[^\n<]*)"
    r"|(?P<STAGE>\[Stage\s*\d+/\d+\][^\n<]*)"
    r"|(?P<SAMPLE>\[(?!(?:INFO|DEBUG|WARNING|ERROR|CRITICAL|Stage|BLOCK|W::|E::|M::|fast-screen|adna|bwa|samtools)\b)[A-Za-z0-9_\-\.]{2,}\])"
    r"|(?P<METRIC>(?:(?<!\w)-?\d+(?:\.\d+)?(?:e[+-]?\d+)?(?:x|%|s|bp|kb|MB|GB)(?!\w))|(?<!\w)-?\d+(?:\.\d+)?e[+-]?\d+\b|(?:\bMAPQ\s*=\s*\d+\b)|(?:\bQ\d{2}\b))"
    r"|(?P<ICON>[✅⚠️🎉🔒❌])"
)

TRACEBACK_START_RE = re.compile(r"^Traceback \(most recent call last\):")
TRACEBACK_FILE_RE = re.compile(r'^\s*File ".*?", line \d+, in .*')
EXCEPTION_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*(?:Error|Exception|Warning):\s*.*")


class LogHighlighter:
    """Parses and formats raw terminal execution logs into safe, styled HTML."""

    @staticmethod
    def strip_ansi(text: str) -> str:
        return ANSI_ESCAPE_RE.sub("", text)

    @staticmethod
    def normalize_carriage_returns(text: str) -> List[str]:
        """Resolves progress bar \r overwrites, isolating the final line state."""
        resolved_lines = []
        for raw_line in text.split("\n"):
            cleaned = LogHighlighter.strip_ansi(raw_line)
            if "\r" in cleaned:
                segments = [seg.strip() for seg in cleaned.split("\r") if seg.strip()]
                if segments:
                    resolved_lines.append(segments[-1])
            else:
                resolved_lines.append(cleaned)
        return resolved_lines

    @classmethod
    def tail_lines_binary(cls, filepath: Path, n_lines: int = 60, chunk_size: int = 4096) -> List[str]:
        """
        Reads backwards from EOF in chunks of 4 KB until n_lines newlines are found.
        Avoids loading large files into memory and eliminates full-file regex operations.
        """
        if not filepath.is_file() or filepath.stat().st_size == 0:
            return []

        try:
            with open(filepath, "rb") as fp:
                fp.seek(0, os.SEEK_END)
                file_size = fp.tell()
                buffer = bytearray()
                lines_found = 0
                pos = file_size

                while pos > 0 and lines_found <= (n_lines + 5):
                    read_len = min(chunk_size, pos)
                    pos -= read_len
                    fp.seek(pos)
                    chunk = fp.read(read_len)
                    buffer = chunk + buffer
                    lines_found = buffer.count(b"\n")

                decoded = buffer.decode("utf-8", errors="replace")
                normalized = cls.normalize_carriage_returns(decoded)
                return normalized[-n_lines:]
        except OSError:
            return []

    @classmethod
    def format_line(cls, line: str, in_traceback: bool = False) -> Tuple[str, bool]:
        """Tokenizes a single line in a single pass while tracking Traceback state."""
        escaped_line = html.escape(line)
        new_in_traceback = in_traceback

        if TRACEBACK_START_RE.match(line):
            new_in_traceback = True
            return f'<span class="ln-error">{escaped_line}</span>', new_in_traceback

        if in_traceback:
            if TRACEBACK_FILE_RE.match(line):
                return f'<span class="ln-traceback">{escaped_line}</span>', True
            if EXCEPTION_NAME_RE.match(line):
                new_in_traceback = False
                return f'<span class="ln-error">{escaped_line}</span>', False
            if line.startswith("    "):
                return f'<span class="ln-traceback">{escaped_line}</span>', True

        is_critical = any(tag in line for tag in ("[CRITICAL]", "CRITICAL:"))
        is_error = any(tag in line for tag in ("[ERROR]", "ERROR:"))

        def _replace_token(match: re.Match) -> str:
            group_dict = match.groupdict()
            for key, val in group_dict.items():
                if val is not None:
                    if key == "TIMESTAMP":
                        return f'<span class="tk-ts">{val}</span>'
                    elif key == "LEVEL":
                        lvl = val.strip("[]")
                        return f'<span class="tk-lvl tk-{lvl.lower()}">[{lvl}]</span>'
                    elif key == "BLOCK":
                        return f'<span class="tk-block">{val}</span>'
                    elif key == "STAGE":
                        return f'<span class="tk-stage">{val}</span>'
                    elif key == "SAMPLE":
                        return f'<span class="tk-sample">{val}</span>'
                    elif key == "METRIC":
                        return f'<span class="tk-metric">{val}</span>'
                    elif key == "ICON":
                        return f'<span class="tk-icon">{val}</span>'
            return match.group(0)

        tokenized = TOKEN_RE.sub(_replace_token, escaped_line)
        if is_critical:
            return f'<span class="ln-critical">{tokenized}</span>', False
        if is_error:
            return f'<span class="ln-error">{tokenized}</span>', False

        return tokenized, new_in_traceback

    @classmethod
    def render_terminal_html(cls, lines: List[str], theme: Optional[Dict[str, str]] = None) -> str:
        """Wraps tokens into a scoped pre block styled according to the active theme palette."""
        formatted_rows = []
        tb_state = False
        for line in lines:
            rendered, tb_state = cls.format_line(line, in_traceback=tb_state)
            formatted_rows.append(rendered)

        body = "\n".join(formatted_rows)

        # Fallback values if theme is None or keys are missing
        code_bg = theme.get("code_bg", "#04090F") if theme else "#04090F"
        code_text = theme.get("code_text", "#FFFFFF") if theme else "#FFFFFF"
        accent = theme.get("accent", "#06B6D4") if theme else "#06B6D4"

        return f"""
        <style>
        .hud-terminal-container {{
            background-color: {code_bg} !important;
            border: 1px solid {accent}55 !important;
            border-radius: 8px;
            padding: 14px 18px;
            font-family: 'Consolas', 'Fira Code', monospace;
            font-size: 0.86rem;
            line-height: 1.5;
            max-height: 540px;
            overflow-y: auto;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.4);
        }}
        .hud-terminal {{
            margin: 0;
            white-space: pre-wrap;
            word-break: break-all;
            color: {code_text} !important;
        }}
        .hud-terminal span {{ font-weight: normal; }}
        .hud-terminal .tk-ts {{ color: #94A3B8 !important; }}
        .hud-terminal .tk-lvl.tk-info {{ color: #38BDF8 !important; font-weight: bold !important; }}
        .hud-terminal .tk-lvl.tk-debug {{ color: #64748B !important; font-weight: bold !important; }}
        .hud-terminal .tk-lvl.tk-warning {{ color: #FBBF24 !important; font-weight: bold !important; }}
        .hud-terminal .tk-lvl.tk-error {{ color: #F87171 !important; font-weight: bold !important; }}
        .hud-terminal .tk-lvl.tk-critical {{
            color: #FFFFFF !important;
            background-color: #DC2626 !important;
            font-weight: 800 !important;
            border-radius: 3px;
            padding: 1px 5px;
            display: inline-block;
            box-shadow: 0 0 6px rgba(220, 38, 38, 0.6);
        }}
        .hud-terminal .tk-sample {{ color: #F472B6 !important; font-weight: bold !important; }}
        .hud-terminal .tk-block {{ color: #34D399 !important; font-weight: bold !important; }}
        .hud-terminal .tk-stage {{ color: #2DD4BF !important; font-weight: bold !important; }}
        .hud-terminal .tk-metric {{ color: #FDE047 !important; font-weight: bold !important; }}
        .hud-terminal .ln-error {{ background: rgba(239, 68, 68, 0.22); color: #FCA5A5 !important; display: block; border-radius: 2px; padding: 0 2px; }}
        .hud-terminal .ln-critical {{
            background: rgba(220, 38, 38, 0.28) !important;
            color: #FEE2E2 !important;
            border-left: 4px solid #EF4444 !important;
            display: block;
            border-radius: 2px;
            padding: 2px 6px;
            margin: 2px 0;
            font-weight: 600;
        }}
        .hud-terminal .ln-traceback {{ color: #F87171 !important; font-style: italic; display: block; }}
        </style>
        <div class="hud-terminal-container">
            <pre class="hud-terminal">{body}</pre>
        </div>
        """