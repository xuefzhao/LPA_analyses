LPA Analyses
============

This repository contains WDL workflows and helper scripts used for analyzing genes in segmental duplicates (eg. LPA).

Project aim
-----------
This project aims to establish scripts and workflows to analyze genes within segmental duplications, focusing on assembly, alignment, and extraction steps that support downstream gene-level analyses in duplicated regions.

Contents
--------
- `wdl/` — WDL workflows for assembly, alignment, extraction, and analysis (examples include `AlignAssemblies.wdl`, `AssemblySeqLpaAnalyses.wdl`, `Paraphase.wdl`, and others).
- Utility scripts — small helpers like `extract_spanned_regions.py` and shell wrappers.

Quick start
-----------
1. Install a WDL runner (Cromwell or miniWDL) and Java/Python as required.
2. From this repository root, run a workflow using your runner, for example:

```bash
# Example with Cromwell
cromwell run wdl/AssemblySeqLpaAnalyses.wdl --inputs inputs.json
```

Notes
-----
- Workflows assume appropriate reference files and inputs are provided via the runner inputs JSON or command-line.
- Check individual WDL files in `wdl/` for workflow-specific requirements and descriptions.

Maintainer
----------
For questions or changes, open an issue or contact the repository maintainer.
