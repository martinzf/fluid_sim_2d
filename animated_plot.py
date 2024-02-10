import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use('fast')

# Animation time resolution
FPS = 40
INTERVAL = 1e3 / FPS

def init(x, y, vmin, vmax):
    plt.title('Vorticity simulation')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.xlim(x[0, 0], x[0, -1])
    plt.ylim(y[0, 0], y[-1, 0])
    z = np.empty(np.shape(x))
    z[0, 0] = vmin
    z[0, 1] = vmax
    lims = plt.contourf(x, y, z)
    cbar = plt.colorbar(lims)
    cbar.set_label('$\omega$ (s$^{-1}$)')
    for child in lims.collections:
        child.remove()

def downsampling(w, dt):
    if 1e3 * dt < INTERVAL:
        indices = np.arange(0, w.shape[0], INTERVAL/(1e3*dt)).round(decimals=0).astype(int)
        return w[indices]
    return w

fig = plt.figure('Simulation', figsize=(7, 6))
cont = plt.contourf(np.empty((2, 2)), np.empty((2, 2)), np.empty((2, 2)))
def animate(x, y, w, dt):
    w_d = downsampling(w, dt)
    vmin, vmax = np.min(w_d), np.max(w_d)
    init_func = lambda : init(x, y, vmin, vmax)
    def func(i):
        global cont
        for child in cont.collections:
            child.remove()
        cont = plt.contourf(x, y, w_d[i], vmin=vmin, vmax=vmax)
    ani = animation.FuncAnimation(fig, func, init_func=init_func, 
                                  frames=int(w_d.shape[0]), interval=INTERVAL, repeat=False)
    #ani.save('preview.gif', fps=FPS)
    plt.show()