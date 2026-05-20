# Bioinformatics Tutorials

This site collects practical, reproducible tutorials for structural bioinformatics workflows. It is organized as a technical handbook: each page is focused on a specific task, toolchain, or execution environment.

The documentation is designed for GitHub Pages and uses a Read the Docs style layout with a persistent sidebar, searchable pages, and compact navigation.

## Tutorial Areas

| Area | Use it for |
| --- | --- |
| Enzyme design | Enzeptional and GT4SD workflows for catalytic feasibility and activity-oriented predictions. |
| Docking | Hydrated docking workflows with explicit waters, AutoDock Vina, Meeko, and AutoGrid4. |
| Molecular dynamics | AMBER and GROMACS protocols for protein-ligand systems, trajectory preparation, and analysis. |
| QM/MM | AMBER + QUICK and GROMACS + CP2K workflows for hybrid quantum/classical simulations. |
| HPC | SLURM job submission, monitoring, and practical execution patterns for simulation workflows. |
| Command-line references | Docker and Git command references for day-to-day project work. |

## How the Documentation Is Organized

The site source lives in `docs/`.

Tutorial pages are stored in `docs/tutorials/` and registered in `mkdocs.yml`. When a new tutorial README is added to the project, translate or write it in English, place the documentation page under `docs/tutorials/`, and add it to the navigation.

## Local Preview

Install the documentation dependencies and start the local server:

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements-docs.txt
mkdocs serve
```

Build the site before publishing:

```bash
mkdocs build --strict
```
