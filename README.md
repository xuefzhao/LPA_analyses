LPA Analyses
============

This repository contains WDL workflows and helper scripts used for analyzing genes in segmental duplicates (eg. LPA).

Project aim
-----------
This project aims to establish scripts and workflows to analyze genes within segmental duplications, focusing on assembly, alignment, and extraction steps that support downstream gene-level analyses in duplicated regions.

Repository layout
-----------------
- `wdl/` — WDL workflow files and workflow-only helper WDLs (`Structs.wdl`, `Utils.wdl`, etc.).
- `scripts/` — Utility scripts used by the workflows (moved here): `extract_spanned_regions.py`, `run_blast_from_table.sh`.
- `examples/` — Example inputs and templates (contains `inputs.template.json`).

Contents
--------
- `wdl/` — WDL workflows for assembly, alignment, extraction, and analysis (examples include `AlignAssemblies.wdl`, `AssemblySeqLpaAnalyses.wdl`, `Paraphase.wdl`, and others).
- Utility scripts — small helpers like `extract_spanned_regions.py` and shell wrappers (now under `scripts/`).


Quick start
-----------
1. Install a WDL runner (Cromwell or miniWDL) and Java/Python as required.
2. From this repository root, run a workflow using your runner, for example:

```bash
# Example with Cromwell
cromwell run wdl/AssemblySeqLpaAnalyses.wdl --inputs examples/inputs.template.json
```

Notes
-----
- Workflows assume appropriate reference files and inputs are provided via the runner inputs JSON or command-line.
- Check individual WDL files in `wdl/` for workflow-specific requirements and descriptions.

Examples and templates
----------------------
See `examples/inputs.template.json` for a minimal example of required inputs and how to point workflows to the helper scripts in `scripts/`.

Maintainer
----------
For questions or changes, open an issue or contact the repository maintainer.
