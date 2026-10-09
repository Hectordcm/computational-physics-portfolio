# 2D and 3D Mass Distributions and Inertia Tensors

This project loads a CSV table of point masses and their coordinates, calculates the
center of mass, constructs the inertia tensor in two and three dimensions, and uses
eigenvalues and eigenvectors to identify principal axes. It produces plots of the
mass distribution and principal directions.

## Methods and tools

- Python, NumPy, Pandas, Matplotlib
- Center-of-mass calculations
- 2D and 3D inertia tensors
- Eigenvalue/eigenvector decomposition with `numpy.linalg.eigh`
- 2D and 3D visualization

## Run

The file `datosmasas.csv` must remain in the same folder as the Python script.

```bash
python -m pip install -r requirements.txt
python projects/inertia-tensor-analysis/inertia_tensor_analysis.py
```

## Data note

The CSV included here is the dataset supplied with this project. Before making the
repository public, confirm that you are allowed to redistribute it (for example, if it
was provided as part of a course assignment). If not, remove the CSV and document how
users can supply their own compatible dataset.
