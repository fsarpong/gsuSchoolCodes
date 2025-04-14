#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 05 2023

@author: felixsarpong
"""

import numpy as np
import matplotlib.pyplot as plt

def solve_ftcs(dx, dt, t_max, order):
    N = int(1 / dx) + 1
    x = np.linspace(0, 1, N)
    r = dt / (dx ** 2)
    U = np.cos(np.pi * x / 2)
    U[-1] = 0
    U_current = U.copy()
    for _ in range(int(t_max / dt)):
        U_new = U_current.copy()
        for j in range(1, N-1):
            U_new[j] = U_current[j] + r * (U_current[j+1] - 2 * U_current[j] + U_current[j-1])
        if order == 1:
            U_new[0] = U_new[1]
        else:
            U_new[0] = U_current[0] + 2 * r * (U_current[1] - U_current[0])
        U_new[-1] = 0
        U_current = U_new
    return x, U_current

# Parameters
params = [(0.1, 0.004), (0.05, 0.001)]
results = {}

for dx, dt in params:
    for order in [1, 2]:
        x, U_num = solve_ftcs(dx, dt, 1, order)
        U_exact = np.exp(-(np.pi / 2) ** 2) * np.cos(np.pi * x / 2)
        results[(dx, order)] = (x, U_num, U_exact)

# Plotting
for dx, dt in params:
    plt.figure()
    for order in [1, 2]:
        x, U_num, U_exact = results[(dx, order)]
        plt.plot(x, U_num, label=f'Order {order}, dx={dx}')
    plt.plot(x, U_exact, '--', label='Exact')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('U')
    plt.savefig(f'problem_II_dx_{dx}.png')
    plt.close()

# Errors
for dx, dt in params:
    for order in [1, 2]:
        _, U_num, U_exact = results[(dx, order)]
        error = np.max(np.abs(U_num - U_exact))
        print(f'dx={dx}, dt={dt}, Order {order} Error: {error}')