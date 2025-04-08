import numpy as np
import matplotlib.pyplot as plt

L = 1
T = 0.1
h = 0.1
k = 0.001
n = int(L / h)
m = int(T / k)
x = np.linspace(0, L, n + 1)
λ = k / h**2

def eksplisitt():
    u = np.zeros((m+1, n+1))
    u[0,:] = np.sin(x)
    u[0,0] = 0
    u[0,-1] = 0
    for j in range(0, m):
        u[j+1, 1:n] = u[j, 1:n] + λ * (u[j, 2:] - 2 * u[j, 1:n] + u[j, :-2])
    return u

def implisitt():
    u = np.zeros((m+1, n+1))
    u[0,:] = np.sin(x)
    u[0,0] = 0
    u[0,-1] = 0
    main = (1 + 2 * λ) * np.ones(n-1)
    off = -λ * np.ones(n-2)
    A = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    for j in range(m):
        b = u[j, 1:n]
        u[j+1, 1:n] = np.linalg.solve(A, b)
    return u

def crank_nicolson():
    u = np.zeros((m+1, n+1))
    u[0,:] = np.sin(x)
    u[0,0] = 0
    u[0,-1] = 0
    main_A = (1 + λ) * np.ones(n-1)
    naboer_A = (-λ/2) * np.ones(n-2)
    A = np.diag(main_A) + np.diag(naboer_A, 1) + np.diag(naboer_A, -1)
    main_B = (1 - λ) * np.ones(n-1)
    naboer_B = (λ/2) * np.ones(n-2)
    B = np.diag(main_B) + np.diag(naboer_B, 1) + np.diag(naboer_B, -1)
    for j in range(0, m):
        b = B @ u[j, 1:n]
        u[j+1, 1:n] = np.linalg.solve(A, b)
    return u

def analytisk_løsning(x, t, N=100):
    u = np.zeros_like(x)
    for n in range(1, N+1, 2):
        bn = (2 * n * np.pi * (-1)**((n+1))) / (n**2 * np.pi**2 - 1)
        u += bn * np.sin(n * np.pi * x) * np.exp(-n**2 * np.pi**2 * t)
    return u

u_eks = eksplisitt()
u_imp = implisitt()
u_cn = crank_nicolson()
u_ana = analytisk_løsning(x, T)

plt.figure(figsize=(10, 5))
plt.plot(x, u_eks[-1], '--', label="Eksplisitt")
plt.plot(x, u_imp[-1], '-.', label="Implisitt")
plt.plot(x, u_cn[-1], ':', label="Crank-Nicolson")
plt.plot(x, u_ana, '-', label="Analytisk", color='black')
plt.xlabel("x")
plt.ylabel("temperatur")
plt.grid(True)
plt.legend()
plt.show()