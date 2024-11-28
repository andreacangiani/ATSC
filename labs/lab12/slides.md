<!--
title: Lab 12
paginate: true

_class: titlepage
-->

# Lab 12
<br>

## Stationary Navier-Stokes equations. deal.II step-57.
<br>

#### Advanced Topic in Scientific Computing - SISSA, UniTS, 2024-2025

###### Pasquale Claudio Africa

###### 28 Nov 2024

---

# Assignment

Read the documentation of `step-57`.

1. It is easy to compare the currently implemented linear solver to just using UMFPACK for the whole linear system. You need to remove the nullspace containing the constant pressures and it is done in `step-56`.
2. For larger computations, especially in 3D, it is necessary to implement MPI parallel solvers and preconditioners. A good starting point would be `step-55`, which uses algebraic multigrid for the velocity block for the Stokes equations. Another option would be to take a look at the list of codes in the [deal.II code gallery](https://dealii.org/developer/doxygen/deal.II/code_gallery_time_dependent_navier_stokes.html), which already contains parallel Navier-Stokes solvers.
