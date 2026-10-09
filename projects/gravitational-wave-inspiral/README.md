# Gravitational-Wave Inspiral Analysis

A scientific-computing project exploring a gravitational-wave strain waveform from the
public Simulating eXtreme Spacetimes (SXS) catalog. The analysis examines the phase,
compares interpolation methods, estimates derivatives, locates an ISCO-related crossing,
and compares numerical solutions of a TaylorT4-inspired evolution equation.

## Methods and tools

- Python, NumPy, Matplotlib, SciPy
- `sxs` package for loading the SXS waveform
- Lagrange interpolation, cubic splines, and PCHIP
- Numerical differentiation and root finding (bisection and Newton–Raphson)
- ODE integration with RK4, RK45, and LSODA
- Error metrics and a convergence comparison for RK4

## Run

From the repository root, install dependencies:

```bash
python -m pip install -r requirements.txt
python projects/gravitational-wave-inspiral/gravitational_wave_inspiral.py
```

The script loads the SXS dataset identified as `SXS:BHNS:0001`. Initial execution may
require an internet connection and a download through the `sxs` package.

## Notes

This is an academic computational project, not a production waveform model. Results depend
on the selected waveform, sampling, interpolation choices, and numerical tolerances.
