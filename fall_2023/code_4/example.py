#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 2023

@author: felixsarpong
"""

import numpy as np
import numerical_solutions as ns
# Parameters
a = 1.0
nu = 0.1
dx = 0.1
dt = 0.005
N = int(1 / dx) + 1
T = 0.1
x = np.linspace(0, 1, N)
U0 = np.sin(2 * np.pi * x)

# Problem 2
is_stable_ftcs = ns.ftcs_stability(a, nu, dx, dt, N)
print(f"FTCS scheme is stable: {is_stable_ftcs}")

# Problem 3
U_lag = ns.btcs_burgers_lag(a, nu, dx, dt, N, T, U0)
U_linear = ns.btcs_burgers_linear(a, nu, dx, dt, N, T, U0)
U_newton = ns.btcs_burgers_newton(a, nu, dx, dt, N, T, U0)
print(f"Burgers (Lag) at t={T}: {U_lag}")
print(f"Burgers (Linear) at t={T}: {U_linear}")
print(f"Burgers (Newton) at t={T}: {U_newton}")

# Problem 4
is_stable_cn = ns.crank_nicolson_stability(nu, dx, dt, N)
print(f"Crank-Nicolson scheme is stable: {is_stable_cn}")