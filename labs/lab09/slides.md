<!--
title: Lab 09
paginate: true

_class: titlepage
-->

# Lab 09
<br>

## Time-dependent problems. deal.II step-26.
<br>

#### Advanced Topic in Scientific Computing - SISSA, UniTS, 2024-2025

###### Pasquale Claudio Africa

###### 14 Nov 2024

---

# Assignment

- Read the documentation of `step-26`.
- Modify `step-26` so to disable local refinement and to utilize explicit Runge-Kutta methods (among the ones available in deal.II) for time integration. Use `step-52` as a reference, specifically examining its `explicit_method`.
- Using the method of manufactured solutions, design an exact solution that satisfies the equations in `step-26`. Compute and visualize the $L^2$ and $H^1$ error w.r.to the exact solution at the final timestep for different types of time discretizations.
