import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

L = 1  #lengde på stang
T = 0.1  #total tid
h = 0.1 #må ha minst tre punkter
k = 0.001

def crank_nicolson():
    n = int(L/h)
    m = int(T/k)
    x = np.linspace(0, L, n+1)
    u = np.zeros((m+1, n+1))
    u[0,:] = np.sin(x)
    u[0,0] = 0
    u[0,-1] = 0

    λ = k / h**2

    #matrise A
    main_A = (1 + λ) * np.ones(n-1)
    naboer_A = (-λ/2) * np.ones(n-2)
    A = np.diag(main_A) + np.diag(naboer_A, 1) + np.diag(naboer_A, -1)

    #matrise B
    main_B = (1 - λ) * np.ones(n-1)
    naboer_B = (λ/2) * np.ones(n-2)
    B = np.diag(main_B) + np.diag(naboer_B, 1) + np.diag(naboer_B, -1)

    for j in range(0, m):
        b = B @ u[j, 1:n] 
        u[j+1, 1:n] = np.linalg.solve(A, b)

    return x, u, m, n

def plott(x, u, m):
    fig, ax = plt.subplots(figsize=(8, 4))
    line, = ax.plot(x, u[0], color="m")
    ax.set_ylim(-1, 1)
    ax.set_xlabel("x")
    ax.set_ylabel("temperatur")

    def update(frame):
        line.set_ydata(u[frame])
        ax.set_title(f"(t = {frame * k:.4f}s)")
        return line,

    ani = animation.FuncAnimation(fig, update, frames=range(0, m+1, max(1, m//200)), interval=20)
    plt.show()

x, u, m, n = crank_nicolson()
plott(x, u, m)