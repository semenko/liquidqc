#!/usr/bin/env python3
"""Capture RustQC CLI output as an SVG using Rich.

Runs RustQC with a deliberate strandedness mismatch to trigger the warning,
then renders the terminal output as an SVG file for documentation.

Usage:
    uv run --with rich docs/scripts/capture_cli_svg.py
"""

import subprocess
import sys

from rich.console import Console
from rich.text import Text

# Paths — adjust these if your test data is elsewhere
RUSTQC = "./target/release/rustqc"
GTF = "../RustQC-benchmarks/test-data/rna/small/chr6.gtf.gz"
BAM = "../RustQC-benchmarks/test-data/rna/small/test.bam"
OUTDIR = "/tmp/rustqc_svg_capture"
SVG_PATH = "docs/public/examples/cli-strandedness-warning.svg"

TITLE = "rustqc rna"


def main():
    # Run RustQC with wrong strandedness to trigger the warning
    cmd = [
        RUSTQC,
        "rna",
        "--stranded",
        "forward",
        "-p",
        "--gtf",
        GTF,
        BAM,
        "--outdir",
        OUTDIR,
    ]
    env = {
        "CLICOLOR_FORCE": "1",
        "PATH": "/usr/bin:/bin:/usr/local/bin",
    }
    result = subprocess.run(
        cmd,
        capture_output=True,
        env=env,
        timeout=120,
    )

    # Combine stdout + stderr (RustQC writes UI to stderr)
    output = (result.stdout + result.stderr).decode("utf-8", errors="replace")

    # Trim to only show the warning box and lines after it.
    # The warning box starts with a line containing "╭─ Warning".
    lines = output.splitlines(keepends=True)
    start = 0
    for i, line in enumerate(lines):
        if "Warning" in line and "╭" in line:
            # Include one blank line before the box for spacing
            start = max(0, i - 1)
            break
    output = "".join(lines[start:]).rstrip() + "\n"

    # Create a Rich console that records output
    console = Console(record=True, width=100)

    # Feed the ANSI output through Rich so it interprets the escape codes
    text = Text.from_ansi(output)
    console.print(text)

    # Export as SVG
    svg = console.export_svg(title=TITLE)
    with open(SVG_PATH, "w") as f:
        f.write(svg)
    print(f"SVG written to {SVG_PATH}")


if __name__ == "__main__":
    main()
