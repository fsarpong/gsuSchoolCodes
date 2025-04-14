# Advanced Numerical Analysis - Homework 4

This repository contains Python implementations of finite difference methods for solving partial differential equations (PDEs) as part of the Advanced Numerical Analysis course (MATH 8610, Fall 2023). This project focuses on:

- Stability analysis of the FTCS scheme for the heat equation (Problem 2),
- Various implementations of the BTCS scheme for the viscous Burgers' equation (Problem 3),
- and stability analysis of the Crank-Nicolson scheme using Gerschgorin's circle theorem (Problem 4).

*Note: Problem 1 has been excluded.*

---

## Problem Descriptions

### Problem 2: Stability Analysis of the FTCS Scheme for the Heat Equation

We consider the following initial-boundary-value problem for the heat equation:
\[
\frac{\partial U}{\partial t} + a \frac{\partial U}{\partial x} = \nu \frac{\partial^2 U}{\partial x^2}, \quad x \in (0,1), \; t > 0
\]
with Dirichlet boundary conditions:
\[
U(0, t) = U(1, t) = 0, \quad t \ge 0.
\]

A matrix-eigenvalue approach is used to investigate the stability of the FTCS (Forward Time Centered Space) scheme.

---

### Problem 3: BTCS Schemes for the Viscous Burgers' Equation

In this problem, we solve the viscous Burgers' equation using a Backward Time Centered Space (BTCS) approach. Three methods are implemented for handling the nonlinear term:
- **Lag Method:** Uses the previous time level to approximate the nonlinear term.
- **Linearization:** Linearizes the nonlinear term about the previous time step.
- **Newton’s Method:** Iteratively solves the nonlinear system using Newton's method.

These are implemented in the functions:
- `btcs_burgers_lag`
- `btcs_burgers_linear`
- `btcs_burgers_newton`

---

### Problem 4: Stability Analysis of the Crank-Nicolson Scheme with Gerschgorin's Circle Theorem

This problem investigates the stability of the Crank-Nicolson scheme for the heat equation under two types of boundary conditions:

1. **Neumann Boundary Conditions:**  
   \[
   U_x(0, t) = a, \quad U_x(1, t) = b, \quad t \geq 0.
   \]
2. **Dirichlet Boundary Conditions:**  
   \[
   U(0, t) = U(1, t) = 0, \quad t \geq 0.
   \]

The stability analysis leverages Gerschgorin's circle theorem to determine bounds on the eigenvalues of the discretized system.

---

## Code Structure

- **example.py**  
  Main driver script that demonstrates the execution of the solvers corresponding to Problems 2, 3, and 4. It prints stability results, convergence rates, and saves plots for visual inspection.

- **numerical_solutions.py**  
  Contains the implementation of the finite difference schemes, including:
  - The FTCS scheme for the heat equation (Problem 2),
  - BTCS variants (lag, linearized, and Newton’s methods) for the Burgers’ equation (Problem 3),
  - and the Crank-Nicolson scheme with the associated stability analysis using Gerschgorin's circle theorem (Problem 4).

- **Homework4.pdf**  
  The assignment document outlining the problems and requirements for Homework 4.

---

## Usage Instructions

### Requirements

- Python 3.x
- NumPy
- SciPy
- Matplotlib

### Running the Code

1. **Clone or download the repository** and navigate to the project directory in your terminal.
2. **Run the main example script:**
   ```bash
   python3 example.py
