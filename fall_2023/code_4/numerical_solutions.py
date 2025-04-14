#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 2023

@author: felixsarpong
"""

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Problem 2: Matrix-Eigenvalue Approach for FTCS Scheme Stability
def ftcs_stability(a, nu, dx, dt, N):
    r = nu * dt / dx**2
    c = a * dt / dx
    A = np.zeros((N-1, N-1))
    for i in range(N-1):
        A[i, i] = 1 - 2 * r
        if i > 0:
            A[i, i-1] = r - c / 2
        if i < N-2:
            A[i, i+1] = r + c / 2
    eigenvalues = np.linalg.eigvals(A)
    return np.max(np.abs(eigenvalues)) <= 1

# Problem 3: BTCS Scheme for Viscous Burgers' Equation
# (i) Lag Nonlinear Term
def btcs_burgers_lag(a, nu, dx, dt, N, T, U0):
    r = nu * dt / dx**2
    U = U0.copy()
    for _ in range(int(T / dt)):
        A = diags([-r, 1 + 2*r, -r], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
        b = U[1:-1] - (dt / (2 * dx)) * U[1:-1] * (U[2:] - U[:-2])
        U[1:-1] = spsolve(A, b)
    return U

# (ii) Linearize About Previous Timestep
def btcs_burgers_linear(a, nu, dx, dt, N, T, U0):
    r = nu * dt / dx**2
    U = U0.copy()
    for _ in range(int(T / dt)):
        A = diags([-r, 1 + 2*r, -r], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
        # Linearize: U_j^{n+1} * (U_{j+1}^n - U_{j-1}^n)
        b = U[1:-1] - (dt / (2 * dx)) * U[1:-1] * (U[2:] - U[:-2])
        U[1:-1] = spsolve(A, b)
    return U

# (iii) Newton's Method
def btcs_burgers_newton(a, nu, dx, dt, N, T, U0, max_iter=10, tol=1e-6):
    r = nu * dt / dx**2
    U = U0.copy()
    for _ in range(int(T / dt)):
        U_new = U.copy()
        for _ in range(max_iter):
            F = U_new[1:-1] + (dt / (2 * dx)) * U_new[1:-1] * (U_new[2:] - U_new[:-2]) - r * (U_new[2:] - 2 * U_new[1:-1] + U_new[:-2]) - U[1:-1]
            if np.max(np.abs(F)) < tol:
                break
            J = diags([
                -r - (dt / (2 * dx)) * U_new[1:-1],
                1 + 2 * r + (dt / (2 * dx)) * (U_new[2:] - U_new[:-2]),
                -r + (dt / (2 * dx)) * U_new[1:-1]
            ], [-1, 0, 1], shape=(N-2, N-2)).tocsc()
            dU = spsolve(J, -F)
            U_new[1:-1] += dU
        U = U_new.copy()
    return U

# Problem 4: Stability of Crank-Nicolson Scheme using Gerschgorin's Circle Theorem
def crank_nicolson_stability(nu, dx, dt, N):
    r = nu * dt / (2 * dx**2)
    A = np.zeros((N-1, N-1))
    B = np.zeros((N-1, N-1))
    for i in range(N-1):
        A[i, i] = 1 + 2 * r
        B[i, i] = 1 - 2 * r
        if i > 0:
            A[i, i-1] = -r
            B[i, i-1] = r
        if i < N-2:
            A[i, i+1] = -r
            B[i, i+1] = r
    # Gerschgorin's theorem on A^-1 B
    for i in range(N-1):
        center = B[i, i] / A[i, i]
        radius = sum(abs(B[i, j] / A[i, i]) for j in range(N-1) if j != i)
        if abs(center) + radius > 1:
            return False
    return True