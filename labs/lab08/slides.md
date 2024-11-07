<!--
title: Lab 08
paginate: true

_class: titlepage
-->

# Lab 08
<br>

## Elasticity equations. deal.II step-8.
<br>

#### Advanced Topic in Scientific Computing - SISSA, UniTS, 2024-2025

###### Pasquale Claudio Africa

###### 08 Nov 2024

---

# Assignment

- Read the documentation of `step-8`.
- Modify `step-8` so to implement the following set of boundary conditions:

$$
\left\{
\begin{aligned}
\mathbf{u} & = \mathbf{0} & \quad & \text{on } \{x = -1\}, \\
\mathbf{\sigma}(\mathbf{u}) \mathbf{n} & = \mathbf{g} & \quad &\text{on } \{x = 1\}, \\
\mathbf{\sigma}(\mathbf{u}) \mathbf{n} & = \mathbf{0} & \quad & \text{elsewhere on } \partial\Omega,
\end{aligned}
\right.
$$
where $\mathbf{\sigma}(\mathbf{u}) = C : \varepsilon (\mathbf{u})$ and $\mathbf{g} = [10, 10]^T$.

These conditions model a plate fixed on its left side, free stress conditions on the top and bottom sides, and a normal traction force on its right side.
