# 2D Ising Model: Metropolis Monte Carlo

This project simulates a square-lattice Ising model with periodic boundary conditions.
A Metropolis update proposes single-spin flips and accepts them according to the
Boltzmann criterion. The script tracks energy and magnetization and visualizes the
evolution and final spin configuration.

## Methods and tools

- Python and NumPy
- Metropolis Monte Carlo
- Periodic boundary conditions
- Energy and magnetization observables
- A moving-average heuristic to estimate a burn-in point
- Matplotlib visualizations

## Run

```bash
python -m pip install -r requirements.txt
python projects/ising-model-monte-carlo/ising2d_monte_carlo.py
```

Parameters such as lattice dimensions, temperature, coupling, number of steps, and random
seed are defined near the top of the script.

## Notes

The burn-in estimate is a simple heuristic based on smoothed energy, not a formal
convergence diagnostic. Physical conclusions should be checked against appropriate
finite-size and temperature-dependent analyses.
