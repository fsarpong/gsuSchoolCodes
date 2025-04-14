#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 05 2023

@author: felixsarpong
"""

import numpy as np
import matplotlib.pyplot as plt
from problem_I_ftcs import exact_solution

nu = 1/6
dx = 0.1
dt = 0.01
N = int(1 / dx) + 1
x = np.linspace(0, 1, N)
r = nu * dt / (dx ** 2)

U0 = np.sin(2 * np.pi * x)
U0[0] = 0
U0[-1] = 0

# First step with FTCS
U1 = U0.copy()
for j in range(1, N-1):
    U1[j] = U0[j] + r * (U0[j+1] - 2 * U0[j] + U0[j-1])

# Leapfrog
U_prev = U0.copy()
U_curr = U1.copy()
for n in range(int(1 / dt)):
    U_next = U_prev.copy()
    for j in range(1, N-1):
        U_next[j] = U_prev[j] + 2 * r * (U_curr[j+1] - 2 * U_curr[j] + U_curr[j-1])
    U_prev = U_curr.copy()
    U_curr = U_next.copy()
    if np.max(np.abs(U_curr)) > 10:  # Check for blow-up
        print("Solution unstable at step", n)
        break

plt.plot(x, U_curr, label='t=1 (leapfrog)')
plt.plot(x, exact_solution(x, 1, nu), '--', label='t=1 (exact)')
plt.legend()
plt.savefig('problem_I_leapfrog.png')
plt.close()

print("Leapfrog scheme is unstable for the heat equation.")