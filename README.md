pyrotini
========

[![GitHub Actions Build Status](https://github.com/DanYev/pyrotini/workflows/CI/badge.svg)](https://github.com/DanYev/pyrotini/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/DanYev/pyrotini/branch/main/graph/badge.svg)](https://codecov.io/gh/DanYev/pyrotini/branch/main)

`pyrotini` glues together [PyRosetta](https://www.pyrosetta.org/) protein design with
[Martini](http://cgmartini.nl/) coarse-grained molecular dynamics (CGMD) in
[GROMACS](https://www.gromacs.org/), and adds analysis tools (Dynamic Flexibility Index, DFI) on top.
It is built for running design → coarse-graining → simulation → analysis pipelines on an HPC cluster
(developed for ASU's SOL cluster).

## What it does

1. **Design** (`tutorial/prt_design.py`) — uses PyRosetta to relax a starting structure and run a
   `FastRelax`/`PackRotamersMover` design loop around a target residue, logging mutations and energies
   for each decoy.
2. **Coarse-graining** (`prepare_files`, `martinize_go`) — fetches a Go contact map for the PDB
   structure (`pyrotini/get_go.py`, automated via Selenium against the RCSU Go-map server), builds the
   GROMACS topology, and runs `martinize2` to produce a Martini 3 Go-model representation.
3. **Simulation** (`solvate`, `energy_minimization`, `heatup`, `equilibration`, `md`) — solvates and
   ionizes the system and drives GROMACS through minimization, heat-up, equilibration, and production MD.
4. **Analysis** (`convert_trajectory`, `get_covariance_matrix`, `get_dfis`, `plot_dfi`) — computes
   per-residue covariance and Dynamic Flexibility Index (DFI) profiles from the trajectory and plots them.

All of the above is exposed as plain functions in `pyrotini/pyrotini.py`, importable as:

```python
import pyrotini as prt
```

## Installation

This package depends on PyRosetta and GROMACS, which are not installable from PyPI/conda-forge directly,
so the workflow below assumes an HPC environment (e.g. SOL) with environment modules.

1. **Load cluster modules:**

    ```bash
    module load mamba
    module load gromacs
    ```

2. **Clone the repository and create the environment:**

   ```bash
   git clone https://github.com/DanYev/pyrotini.git
   cd pyrotini
   mamba env create -n prttest --file env.yml
   source activate prttest
   ```

3. **Install the package (editable):**

    ```bash
    pip install -e .
    ```

4. **Install PyRosetta separately** following the instructions at
   [pyrosetta.org](https://www.pyrosetta.org/downloads) — it is not bundled in `env.yml`/`pyproject.toml`
   because it requires an academic/commercial license.

## Usage

The `tutorial/` directory contains end-to-end examples:

- `prt_design.py` — runs a PyRosetta design job around a given residue:

  ```bash
  python tutorial/prt_design.py -f 1btl.pdb -t 44A -r 1
  ```

- `prt_run_all.py` — runs the full CGMD pipeline (file prep → Go-map → Martini coarse-graining →
  solvation → minimization → heat-up → equilibration → production MD → DFI analysis) on a PDB file:

  ```bash
  python tutorial/prt_run_all.py -f 1btl.pdb -d test
  ```

- `submit2sol.sh` / `submitdesign.sh` / `paralleldesign.sh` — SLURM batch scripts for submitting the
  above as jobs on the SOL cluster:

  ```bash
  sbatch tutorial/submit2sol.sh
  ```

Sample input data (`1btl.pdb`, Martini force-field `.itp` files, GROMACS `.mdp` parameter files) ships
under `pyrotini/data/` and is resolved at runtime via `pyrotini.PRT_DICT`.

## Repository layout

```
pyrotini/
├── pyrotini/             # package source
│   ├── pyrotini.py       # design, coarse-graining, MD, and DFI analysis functions
│   ├── get_go.py         # Selenium-based Go contact-map fetcher
│   └── data/             # bundled Martini/GROMACS parameter and structure files
├── tutorial/             # example scripts and SLURM submission scripts
├── devtools/             # conda environments and dev scripts
├── docs/                 # Sphinx documentation
└── env.yml               # conda/mamba environment for non-PyRosetta dependencies
```

## Testing

```bash
pip install -e ".[test]"
pytest pyrotini/tests/
```

## Documentation

Full API docs are built with Sphinx from `docs/` and published to
[pyrotini.readthedocs.io](https://pyrotini.readthedocs.io/).

## Copyright

Copyright (c) 2024, DY

### Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
