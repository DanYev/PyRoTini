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
    module load gromacs
    mamba env create -n prttest --file env.yml
    ```

    activate environment:

    ```bash
    source activate prttest
    ```

2. **Clone this repository:**

   ```bash 
   git clone https://github.com/DanYev/pyrotini.git
   ```

3. **Install the package:**

    ```bash
    pip install -e .
    ```

**Submission:**

```bash
sbatch submit2sol.sh
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
