<!--
title: Lab 10
paginate: true

_class: titlepage
-->

# Lab 10
<br>

## Mixed FEM. deal.II step-20.
<br>

#### Advanced Topic in Scientific Computing - SISSA, UniTS, 2024-2025

###### Pasquale Claudio Africa

###### 21 Nov 2024

---

# Assignment

Read the documentation of `step-20`.

1. Modify `step-20` so to use matrices and vectors from the `TrilinosWrappers` namespace.
2. In the `solve()` method, try different solvers and preconditioners in order to detect the combination minimizing the number of linear solver iterations. Focus in particular on Algebraic Multi-Grid (AMG) preconditioners.
