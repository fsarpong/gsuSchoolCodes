# Heat Equation Numerical Solutions

This project provides numerical solutions to two initial-boundary-value problems for the heat equation using various finite difference methods. The problems are solved using Python with NumPy and Matplotlib for computations and visualizations.

## Problem Description

We consider two problems for the heat equation:

### Problem I

\[
\begin{cases}
\frac{\partial U}{\partial t} = \nu \frac{\partial^{2} U}{\partial x^{2}}, & x \in (0,1),\quad t > 0, \\
U(x, 0) = \sin(2\pi x), & x \in [0,1], \\
U(0, t) = U(1, t) = 0, & t \geq 0.
\end{cases}
\]

### Problem II

\[
\begin{cases}
\frac{\partial U}{\partial t} = \frac{\partial^{2} U}{\partial x^{2}}, & x \in (0,1),\quad t > 0, \\
U(x, 0) = \cos\left(\frac{\pi x}{2}\right), & x \in [0,1], \\
\frac{\partial U}{\partial x}(0, t) = 0, \\
U(1, t) = 0, & t \geq 0.
\end{cases}
\]

## Solution Overview

### (a) Exact Solutions

- **Problem I:** Using separation of variables, the exact solution is:
  
  \[
  U(x, t) = e^{-(2\pi)^2 \nu t} \sin(2\pi x)
  \]
  
  With \(\nu = \frac{1}{6}\), this becomes:
  
  \[
  U(x, t) = e^{-\frac{2\pi^2}{3} t} \sin(2\pi x)
  \]

- **Problem II:** The exact solution is:

  \[
  U(x, t) = e^{-\left(\frac{\pi}{2}\right)^2 t} \cos\left(\frac{\pi x}{2}\right)
  \]

### (b) Second-Order Forward Difference Formula

We derive a second-order accurate forward difference approximation for \( u'(x_j) \):

\[
u'(x_j) \approx \frac{-3 u(x_j) + 4 u(x_j + \Delta x) - u(x_j + 2\Delta x)}{2 \Delta x}
\]

### (c) Numerical Solution for Problem I

- **Method:** Forward Time Centered Space (FTCS) scheme.
- **Parameters:** \(\nu = \frac{1}{6}\), \(\Delta x = 0.1\), \(\Delta t = 0.01\).
- **Results:** Solutions are computed at \( t = 0.01,\, 0.1,\, 1,\, 10 \) and compared with the exact solution. The largest stable \(\Delta t\) is approximately 0.03.

### (d) Leapfrog Scheme for Problem I

- **Method:** Leapfrog scheme.
- **Observation:** The leapfrog scheme is unstable for the heat equation, leading to rapid growth in the solution.

### (e) Numerical Solution for Problem II

- **Methods:** FTCS with first-order and second-order discretizations of the Neumann boundary condition.
- **Parameters:**
  - Case 1: \(\Delta x = 0.1\), \(\Delta t = 0.004\)
  - Case 2: \(\Delta x = 0.05\), \(\Delta t = 0.001\)
- **Results:** Second-order discretizations provide more accurate solutions, with errors decreasing as \(\Delta x\) decreases.

### Additional Problems

- **Problem 2:** Numerical solution for a heat equation with mixed boundary conditions.
  - **Consistency Analysis:**
    - BTCS is pointwise consistent for the heat equation.
    - FTCS with a first-order Neumann condition is inconsistent unless \( u' = 0 \).

## Code Structure

- **`problem_I_ftcs.py`**: Solves Problem I using FTCS and compares with the exact solution.
- **`problem_I_leapfrog.py`**: Attempts to solve Problem I using the leapfrog scheme, demonstrating instability.
- **`problem_II_neumann.py`**: Solves Problem II with first and second-order Neumann discretizations.
- **`problem_2.py`**: Solves an additional problem with mixed boundary conditions.

## Usage

1. **Requirements:**  
   Ensure Python 3.13.3 is installed along with [NumPy](https://numpy.org/) and [Matplotlib](https://matplotlib.org/).

2. **Running the Scripts:**  
   Run the individual Python scripts to generate solutions and plots. For example:
   ```bash
   python3 problem_I_ftcs.py
