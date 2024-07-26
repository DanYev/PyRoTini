pyrotini
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/pyrotini/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/pyrotini/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pyrotini/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pyrotini/branch/main)


## CGMD scripts for pyRosetta on SOL cluster

### Installation

1. **Install the environment:**

    ```bash
    module load mamba
    mamba env create -n prttest --file env.yml
    ```

    make sure you have conda initialized in your environment:

    ```bash
    conda init bash
    ```

2. **Clone this repository:**

   ```bash
   git clone https://github.com/John-Kazan/OpenMMScripts.git
   ```

3. **Run the setup script:**

    ```bash
    ./setup.sh
    ```

### Running Simulations

#### Using `openmm_npt.py`

This script runs the NPT simulation directly using OpenMM. It utilizes GPU (CUDA) for acceleration.

**Usage:**

```bash
python openmm_npt.py -pdb ./pdb/1btl.pdb
```

Replace `pdb/1btl.pdb` with the path to your PDB file.

#### Using `openmm_npt.sh` (for PHX cluster)

This wrapper script automatically loads the correct CUDA version and OpenMM environment on the PHX cluster.

**Usage:**

```bash
openmm_npt.sh -pdb ./pdb/1btl.pdb
```

Replace `1btl.pdb` with your PDB file or use the one provided for testing.

### Submitting to Scheduler

To submit the simulation to a scheduler (e.g., SLURM), use the `submit_sbatch.sh` script. 

**Before submission:**

Edit `submit_sbatch.sh`:
- Modify the `openmm_npt.sh -pdb ./pdb/1btl.pdb` line to reflect your own PDB.

**Submission:**

```bash
sbatch submit_sbatch.sh
```

To continue the simulation from the last point simply submit the job again. It will automatially continue.

### Notes

- The simulation is set to run for 100ns. On RTX 2080 the performance is ~200ns/day.
- This repository is designed for running simulations on the PHX cluster, but can be adapted for other environments.
- The provided scripts are examples and may need modifications to suit your specific needs.
- For more detailed information about OpenMM, refer to the official documentation: [https://openmm.org/](https://openmm.org/)

### Copyright

Copyright (c) 2024, DY


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
