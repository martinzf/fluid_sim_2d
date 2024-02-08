import numpy as np
import os
import gfort2py as gf 

# Importing Fortran library
SHARED_LIB_NAME = os.path.join(os.getcwd(), 'fluid_sim', 'fluid_sim.dll')
MOD_FILE_NAME = os.path.join(os.getcwd(), 'fluid_sim', 'fluid_sim.mod')
fluid_sim = gf.fFort(SHARED_LIB_NAME, MOD_FILE_NAME)

Nx = 256
Ny = 256
tf = 10
dt = .1
steps = int(tf / dt)
t = np.arange(0, tf, dt)
nu = .1
Lx = 2 * np.pi
Ly = np.pi

x = np.linspace(-Lx/2, Lx/2, Nx)
y = np.linspace(-Ly/2, Ly/2, Ny)
X, Y = np.meshgrid(x, y)

w0 = np.cos(X) * np.cos(Y)
w_cum = np.empty((len(t), *np.shape(w0)))

w = fluid_sim.vort_solve(Nx, Ny, Lx, Ly, nu, dt, steps, w0)[0]
print(w)