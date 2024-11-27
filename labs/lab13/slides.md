<!--
title: Lab 13
paginate: true

_class: titlepage
-->

# Lab 13
<br>

## Time-dependent Navier-Stokes equations. deal.II step-35.
<br>

#### Advanced Topic in Scientific Computing - SISSA, UniTS, 2024-2025

###### Pasquale Claudio Africa

###### 29 Nov 2024

---

# Assignment

Read the documentation of `step-35`.

This program can be extended in the following directions:
1. Adaptive mesh refinement: Using adaptive mesh refinement can lead to increased accuracy while not significantly increasing the computational time.
2. High Reynolds numbers: As we can see from the results, increasing the Reynolds number changes significantly the behavior of the discretization scheme. Using well-known stabilization techniques we could be able to compute the flow in this, or many other problems, when the Reynolds number is very large and where computational costs demand spatial resolutions for which the flow is only marginally resolved, especially for 3D turbulent flows.
