#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 2023

@author: felixsarpong
"""

import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x, t, nu):
    return np.exp(- (2 * np.pi)**2 * nu * t) * np.sin(2 * np.pi * x)

def ftcs_solve(dx, dt, t_max, nu):
    N = int(1 / dx) + 1
    x = np.linspace(0, 1, N)
    r = nu * dt / (dx ** 2)
    U = np.sin(2 * np.pi * x)
    U[0] = 0
    U[-1] = 0
    t = 0
    while t < t_max - dt/2:  # Adjust for floating-point precision
        U_new = U.copy()
        for j in range(1, N-1):
            U_new[j] = U[j] + r * (U[j+1] - 2 * U[j] + U[j-1])
        U = U_new
        t += dt
    return x, U

nu = 0.1
t_max = 0.1
params = [
    (0.1, 0.05, '$\Delta x = 0.1, \Delta t = 0.05$'),
    (0.05, 0.0125, '$\Delta x = 0.05, \Delta t = 0.0125$'),
    (0.025, 0.003125, '$\Delta x = 0.025, \Delta t = 0.003125$')
]

plt.figure(figsize=(10, 6))
errors = []
for dx, dt, label in params:
    x, U_num = ftcs_solve(dx, dt, t_max, nu)
    U_exact = exact_solution(x, t_max, nu)
    plt.plot(x, U_num, label=f'Numerical: {label}')
    error = np.max(np.abs(U_num - U_exact))
    errors.append(error)
    print(f'dx={dx}, dt={dt}, Error={error}')

x_fine = np.linspace(0, 1, 1000)
plt.plot(x_fine, exact_solution(x_fine, t_max, nu), 'k--', label='Exact')
plt.xlabel('x')
plt.ylabel('U(x, 0.1)')
plt.title('FTCS Scheme at t = 0.1')
plt.legend()
plt.grid(True)
plt.savefig('ftcs_solutions.png')
plt.close()

for i in range(1, len(errors)):
    rate_t = np.log(errors[i-1] / errors[i]) / np.log(4)
    rate_x = np.log(errors[i-1] / errors[i]) / np.log(2)
    print(f'Step {i}: Temporal rate = {rate_t:.2f}, Spatial rate = {rate_x:.2f}')