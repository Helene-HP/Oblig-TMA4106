import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


L = 1 #lengde på stang
T = 0.1 #tid
h = 0.1
k = 0.001

def eksplisitt_skjema():
    n = int(L/h)
    m = int(T/k)
    x = np.linspace(0, L, n+1)
    u = np.zeros((m+1, n+1))
    u[0,:] = np.sin(x)
    u[0,0] = 0
    u[0,-1] = 0

    for j in range(0, m):
        u[j+1, 1:n] = u[j, 1:n] + (k / h**2) * (u[j, 2:] - 2 * u[j, 1:n] + u[j, :-2])
    
    return x, u, m, n

def plott(x, u, m):
    fig, ax = plt.subplots(figsize=(8, 4))
    line, = ax.plot(x, u[0], color="m")
    ax.set_ylim(-1, 1)
    ax.set_xlabel("x")
    ax.set_ylabel("temperatur")

    def update(frame):
        line.set_ydata(u[frame])
        ax.set_title(f" (t = {frame * k:.4f}s)")
        return line,

    ani = animation.FuncAnimation(fig, update, frames=range(0, m+1, max(1, m//200)), interval=20)
    plt.show()

x,u,m,n = eksplisitt_skjema()
plott(x,u,m)


    