import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
matplotlib.use('Qt5Agg')
plt.style.use('fast')

# Animation time resolution
FPS = 40

def init(x, y, vmin, vmax):
    plt.xlim(x[0, 0], x[0, -1])
    plt.ylim(y[0, 0], y[-1, 0])
    z = np.empty(np.shape(x))
    z[0, 0] = vmin
    z[0, 1] = vmax
    lims = plt.contourf(x, y, z)
    cbar = plt.colorbar(lims)
    cbar.set_label('$\omega$ (s$^{-1}$)')
    for c in lims.collections:
        c.remove()
    lims = []

def downsampling(w, dt):
    if dt < 1 / FPS:
        indices = np.arange(0, w.shape[0], 1/(FPS*dt)).round(decimals=0).astype(int)
        return w[indices], 1e3 / FPS
    elif dt > 1 / FPS:
        return w, 1e3 * dt
    return w, 1e3 / FPS

fig = plt.figure('Simulation', figsize=(7, 6))
plt.title('Vorticity simulation')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
cont = plt.contourf(np.empty((2, 2)), np.empty((2, 2)), np.empty((2, 2)))
def animate(x, y, w, dt):
    w_d, interval = downsampling(w, dt)
    vmin, vmax = np.min(w_d), np.max(w_d)
    init(x, y, vmin, vmax)
    def func(i):
        global cont
        for c in cont.collections:
            c.remove()
        cont = plt.contourf(x, y, w_d[i], vmin=vmin, vmax=vmax)
        return cont.collections
    ani = animation.FuncAnimation(fig, func, frames=w_d.shape[0], interval=interval, blit=True, repeat=False)
    plt.show()