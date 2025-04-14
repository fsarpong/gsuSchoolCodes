#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 2023

@author: felixsarpong
"""

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def exact_solution(x, t, nu):
    return np.exp(- (2 * np.pi)**2 * nu * t) * np.sin(2 * np.pi * x)

def btcs_solve(dx, dt, t_max, nu):
    N = int(1 / dx) + 1
    x = np.linspace(0, 1, N)
    r = nu * dt / (dx ** 2)
    A = diags([-r, 1 + 2*r, -r], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
    U = np.sin(2 * np.pi * x)
    U[0] = 0
    U[-1] = 0
    t = 0
    while t < t_max - dt/2:
        b = U[1:-1].copy()
        U[1:-1] = spsolve(A, b)
        t += dt
    return x, U

def crank_nicolson_solve(dx, dt, t_max, nu):
    N = int(1 / dx) + 1
    x = np.linspace(0, 1, N)
    r = nu * dt / (2 * dx ** 2)
    A = diags([-r, 1 + 2*r, -r], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
    B = diags([r, 1 - 2*r, r], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
    U = np.sin(2 * np.pi * x)
    U[0] = 0
    U[-1] = 0
    t = 0
    while t < t_max - dt/2:
        b = B @ U[1:-1]
        U[1:-1] = spsolve(A, b)
        t += dt
    return x, U

nu = 0.1
t_max = 0.1
params = [
    (0.1, 0.05, '$\Delta x = 0.1, \Delta t = 0.05$'),
    (0.05, 0.025, '$\Delta x = 0.05, \Delta t = 0.025$'),
    (0.025, 0.0125, '$\Delta x = 0.025, \Delta t = 0.0125$')
]

for scheme, name in [(btcs_solve, 'BTCS'), (crank_nicolson_solve, 'Crank-Nicolson')]:
    plt.figure(figsize=(10, 6))
    errors = []
    for dx, dt, label in params:
        x, U_num = scheme(dx, dt, t_max, nu)
        U_exact = exact_solution(x, t_max, nu)
        plt.plot(x, U_num, label=f'Numerical: {label}')
        error = np.max(np.abs(U_num - U_exact))
        errors.append(error)
        print(f'{name}: dx={dx}, dt={dt}, Error={error}')
    
    x_fine = np.linspace(0, 1, 1000)
    plt.plot(x_fine, exact_solution(x_fine, t_max, nu), 'k--', label='Exact')
    plt.xlabel('x')
    plt.ylabel('U(x, 0.1)')
    plt.title(f'{name} Scheme at t = 0.1')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{name.lower().replace("-", "_")}_solutions.png')
    plt.close()
    
    for i in range(1, len(errors)):
        rate = np.log(errors[i-1] / errors[i]) / np.log(2)
        print(f'{name} Step {i}: Convergence rate = {rate:.2f}')