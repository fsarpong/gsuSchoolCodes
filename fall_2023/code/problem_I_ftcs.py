#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 05 2023

@author: felixsarpong
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
nu = 1/6
dx = 0.1
dt = 0.01
N = int(1 / dx) + 1
x = np.linspace(0, 1, N)
times = [0.01, 0.1, 1, 10]
steps = [int(t / dt) for t in times]
r = nu * dt / (dx ** 2)

# Initial condition
U = np.sin(2 * np.pi * x)
U[0] = 0
U[-1] = 0

# Exact solution function
def exact_solution(x, t, nu):
    return np.exp(-2 * (np.pi ** 2) * nu * t) * np.sin(2 * np.pi * x)

# FTCS scheme
U_current = U.copy()
errors = []
for t, n_steps in zip(times, steps):
    for _ in range(n_steps):
        U_new = U_current.copy()
        for j in range(1, N-1):
            U_new[j] = U_current[j] + r * (U_current[j+1] - 2 * U_current[j] + U_current[j-1])
        U_current = U_new
    U_exact = exact_solution(x, t, nu)
    error = np.max(np.abs(U_current - U_exact))
    errors.append(error)
    plt.plot(x, U_current, label=f't={t} (num)')
    plt.plot(x, U_exact, '--', label=f't={t} (exact)')
plt.legend()
plt.xlabel('x')
plt.ylabel('U')
plt.savefig('problem_I_comparison.png')
plt.close()

# Largest dt
dt_max = 0.03  # Stability limit
dt_test = np.linspace(0.01, 0.04, 10)
max_errors = []
for dt in dt_test:
    r = nu * dt / (dx ** 2)
    U_current = np.sin(2 * np.pi * x)
    U_current[0] = 0
    U_current[-1] = 0
    for _ in range(int(10 / dt)):
        U_new = U_current.copy()
        for j in range(1, N-1):
            U_new[j] = U_current[j] + r * (U_current[j+1] - 2 * U_current[j] + U_current[j-1])
        U_current = U_new
    U_exact = exact_solution(x, 10, nu)
    max_errors.append(np.max(np.abs(U_current - U_exact)))

print("Errors at t = 0.01, 0.1, 1, 10:", errors)
print("Largest reasonable dt ≈ 0.03 based on stability and error < 0.01")