# Computational Physics Portfolio

**Author:** Hector Dario Cortes Martinez  
**Focus:** Computational physics, numerical methods, scientific Python, and quantitative analysis  
**Education:** Physics undergraduate, Universidad Industrial de Santander (UIS), Colombia

This repository collects selected academic projects written in Python. The projects
demonstrate numerical modeling, scientific data analysis, Monte Carlo simulation, and
linear-algebra methods applied to physical systems.

## Selected projects

### 1. Gravitational-Wave Inspiral Analysis
Uses a waveform from the public SXS catalog to explore phase interpolation, numerical
differentiation, ISCO-related crossing estimates, and comparisons between RK4, RK45,
LSODA, and an SXS-derived reference.
- Folder: [`projects/gravitational-wave-inspiral/`](projects/gravitational-wave-inspiral/)

### 2. 2D Ising Model Monte Carlo
Implements Metropolis single-spin updates on a periodic lattice and tracks energy and
magnetization across Monte Carlo steps.
- Folder: [`projects/ising-model-monte-carlo/`](projects/ising-model-monte-carlo/)

### 3. 2D and 3D Inertia Tensor Analysis
Calculates centers of mass and inertia tensors for point-mass distributions, then
obtains principal axes from eigenvectors and visualizes the results.
- Folder: [`projects/inertia-tensor-analysis/`](projects/inertia-tensor-analysis/)

## Requirements

- Python 3.10 or newer recommended
- Dependencies listed in [`requirements.txt`](requirements.txt)

Install them with:

```bash
python -m pip install -r requirements.txt
```

Follow the README inside each project folder for specific run instructions.

## About this portfolio

These are academic computational projects intended to demonstrate scientific reasoning
and programming skills. They are not presented as professional research software.
Some scripts originated in Google Colab and have been adapted for local execution.

## Data and attribution

The gravitational-wave project loads data from the SXS catalog through the `sxs` package.
The inertia-tensor project includes a CSV dataset; confirm permission to redistribute
course-provided data before publishing this repository publicly.
