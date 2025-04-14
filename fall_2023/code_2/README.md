
# Numerical Solutions for PDEs: Heat and Hyperbolic Equations

This project provides numerical solutions to initial-boundary-value problems for hyperbolic and parabolic partial differential equations (PDEs) using finite difference methods. The solutions are implemented in Python with NumPy and Matplotlib for computations and visualizations.

---

## Project Overview

This repository contains Python numnerical implementations for solving two distinct problems related to PDEs:

1. **Explicit Scheme for the Heat Equation**  
   - Solves the heat equation using the FTCS scheme with Dirichlet boundary conditions.  
   - Analyzes temporal and spatial convergence rates.

2. **Implicit Schemes for the Heat Equation**  
   - Implements Backward Time Centered Space (BTCS) and Crank-Nicolson schemes.  
   - Compares their convergence rates with the explicit FTCS scheme.

---

## Problem Descriptions

### Problem 2: Heat Equation with Explicit Scheme
Solve the heat equation:  
\[
\frac{\partial U}{\partial t} = \nu \frac{\partial^2 U}{\partial x^2}, \quad x \in (0,1), \, t > 0
\]
- Initial condition: \( U(x, 0) = \sin(2\pi x) \)  
- Boundary conditions: \( U(0, t) = U(1, t) = 0 \)  
- Diffusion coefficient: \( \nu = 0.1 \)  
- Use the FTCS scheme and compute convergence rates at \( t = 0.1 \) with varying grid sizes.

### Problem 3: Heat Equation with Implicit Schemes
Solve the same heat equation as Problem 2 using:  
- The BTCS scheme.  
- The Crank-Nicolson scheme.  
Compare their convergence rates with the FTCS scheme at \( t = 0.1 \).

---

## Solution Highlights

### Problem 2
- FTCS scheme implemented for the heat equation.  
- Convergence: First-order in time, second-order in space (verified via sup-norm error at \( t = 0.1 \)).

### Problem 3
- BTCS: First-order in time, second-order in space.  
- Crank-Nicolson: Second-order in both time and space.  
- Visual and quantitative comparisons provided.

---

## Code Structure

- **Question2.py**: FTCS scheme for the heat equation.  
  - Outputs: `ftcs_solutions.png` (plots), convergence rates printed.  
- **Question3.py**: BTCS and Crank-Nicolson schemes.  
  - Outputs: `btcs_solutions.png`, `crank_nicolson_solutions.png`, convergence rates printed.

---

## Usage Instructions

### Requirements
- Python 3.x  
- NumPy  
- Matplotlib  
- SciPy (for implicit schemes)

### Running the Code
1. Run `Question2.py` for Problem 2 (FTCS solution and plot).  
2. Run `Question3.py` for Problem 3 (BTCS and Crank-Nicolson solutions and plots).  

### Output
- Plots saved as PNG files.  
- Errors and convergence rates printed to the console.

---

## Key Observations

- **Stability**:  
  - FTCS: Unstable for hyperbolic PDEs; stable for heat equation if \( \frac{\nu \Delta t}{\Delta x^2} \leq \frac{1}{2} \).  
  - BTCS and Crank-Nicolson: Unconditionally stable for the heat equation.  

- **Convergence**:  
  - Heat equation schemes match theoretical convergence rates.

- **Accuracy**:  
  - Crank-Nicolson outperforms BTCS due to second-order temporal accuracy.  
  - FTCS is simpler but less accurate for large time steps.

---

## Conclusion

This project showcases finite difference methods for solving PDEs, emphasizing stability, convergence, and accuracy. The Python code enables replication and experimentation with various parameters and schemes.

