import gui
import animated_plot

# Importing Fortran library
import os
import gfort2py as gf 
SHARED_LIB_NAME = os.path.join(os.getcwd(), 'fluid_sim', 'fluid_sim.dll')
MOD_FILE_NAME = os.path.join(os.getcwd(), 'fluid_sim', 'fluid_sim.mod')
fluid_sim = gf.fFort(SHARED_LIB_NAME, MOD_FILE_NAME)

if __name__ == '__main__':
    Nx, Ny, Lx, Ly, t, dt, steps, nu, x, y, w0 = gui.setup()
    w = fluid_sim.vort_solve(Nx, Ny, Lx, Ly, nu, dt, steps, w0)[0]
    animated_plot.animate(x, y, w, dt)