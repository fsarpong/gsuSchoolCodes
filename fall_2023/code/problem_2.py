#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 05 2023

@author: felixsarpong
"""

import numpy as np
import matplotlib.pyplot as plt

dx = 0.1
dt = 0.004
N = int(1 / dx) + 1
x = np.linspace(0, 1, N)
r = dt / (dx ** 2)

U = np.cos(np.pi * x / 2)
U_current = U.copy()
t = 0
for n in range(int(1 / dt)):
    t = (n + 1) * dt
    U_new = U_current.copy()
    U_new[-1] = np.sin(2 * np.pi * t)
    U_new[0] = U_current[0] + r * (2 * U_current[1] - 2 * U_current[0] - 4 * np.pi * dx)
    for j in range(1, N-1):
        U_new[j] = U_current[j] + r * (U_current[j+1] - 2 * U_current[j] + U_current[j-1])
    U_current = U_new

plt.plot(x, U_current, label='t=1 (num)')
plt.legend()
plt.savefig('problem_2.png')
plt.close()