import numpy as np
import tkinter as tk
import customtkinter as ctk

NVALS = ['256', '512', '1024']
WIDTH, HEIGHT = 550, 300 # pixels

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        # Styling
        self.title('Simulation parameters')    
        self.geometry(f'{WIDTH}x{HEIGHT}') 
        self.frame = ctk.CTkFrame(self) 
        self.frame.pack(padx=20, pady=20, fill='both', expand=True)
        self.frame.columnconfigure(0, weight=1, minsize=100)
        self.frame.columnconfigure(1, weight=1, minsize=100)
        self.frame.columnconfigure(2, weight=1, minsize=100)
        self.frame.columnconfigure(3, weight=1, minsize=100)
        self.frame.columnconfigure(4, weight=1, minsize=100)
        self.frame.columnconfigure(5, weight=1, minsize=100)

        # Nx
        self.NxLabel = ctk.CTkLabel(self.frame, text='Nₓ')
        self.NxLabel.grid(row=0, column=0, padx=10, pady=10, sticky='ew')
        self.NxEntry = tk.StringVar(value=NVALS[0])
        self.NxBox = ctk.CTkComboBox(self.frame, values=NVALS, variable=self.NxEntry)
        self.NxBox.set(NVALS[0])
        self.NxBox.grid(row=0, column=1, columnspan=2, padx=10, pady=10, sticky='ew')

        # Ny
        self.NyLabel = ctk.CTkLabel(self.frame, text='Nᵧ')
        self.NyLabel.grid(row=0, column=3, padx=10, pady=10, sticky='ew')
        self.NyEntry = tk.StringVar(value=NVALS[0])
        self.NyBox = ctk.CTkComboBox(self.frame, values=NVALS, variable=self.NyEntry)
        self.NyBox.set(NVALS[0])
        self.NyBox.grid(row=0, column=4, columnspan=2, padx=10, pady=10, sticky='ew')

        # Lx 
        self.LxLabel = ctk.CTkLabel(self.frame, text='Lₓ (m)')
        self.LxLabel.grid(row=2, column=0, padx=10, pady=10, sticky='ew')
        self.LxEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected positive float', 
                              placeholder_text_color='red')
        self.LxEntry.insert(0, '2*np.pi')
        self.LxEntry.grid(row=2, column=1, columnspan=2, padx=10, pady=10, sticky='ew')

        # Ly 
        self.LyLabel = ctk.CTkLabel(self.frame, text='Lᵧ (m)')
        self.LyLabel.grid(row=2, column=3, padx=10, pady=10, sticky='ew')
        self.LyEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected positive float', 
                              placeholder_text_color='red')
        self.LyEntry.insert(0, '2*np.pi')
        self.LyEntry.grid(row=2, column=4, columnspan=2, padx=10, pady=10, sticky='ew')

        # t 
        self.tLabel = ctk.CTkLabel(self.frame, text='t (s)')
        self.tLabel.grid(row=3, column=0, padx=10, pady=10, sticky='ew')
        self.tEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected positive float', 
                              placeholder_text_color='red')
        self.tEntry.insert(0, '20')
        self.tEntry.grid(row=3, column=1, columnspan=2, padx=10, pady=10, sticky='ew')

        # dt 
        self.dtLabel = ctk.CTkLabel(self.frame, text='dt (s)')
        self.dtLabel.grid(row=3, column=3, padx=10, pady=10, sticky='ew')
        self.dtEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected positive float', 
                              placeholder_text_color='red')
        self.dtEntry.insert(0, '0.01')
        self.dtEntry.grid(row=3, column=4, columnspan=2, padx=10, pady=10, sticky='ew')

        # nu 
        self.nuLabel = ctk.CTkLabel(self.frame, text='ν (m²/s)')
        self.nuLabel.grid(row=4, column=0, padx=10, pady=10, sticky='ew')
        self.nuEntry = ctk.CTkEntry(self.frame, 
                              placeholder_text='Expected float', 
                              placeholder_text_color='red')
        self.nuEntry.insert(0, '1')
        self.nuEntry.grid(row=4, column=1, columnspan=2, padx=10, pady=10, sticky='ew')

        # w0
        self.w0Label = ctk.CTkLabel(self.frame, text='ω₀ (s⁻¹)')
        self.w0Label.grid(row=4, column=3, padx=10, pady=10, sticky='ew')
        self.w0Entry = ctk.CTkEntry(self.frame, 
                                   placeholder_text='Expected real NumPy (Ny, Nx) ndarray', 
                                   placeholder_text_color='red')
        self.w0Entry.insert(0, 'np.cos(x)*np.cos(y)')
        self.w0Entry.grid(row=4, column=4, columnspan=2, padx=10, pady=10, sticky='ew')

        # Simulate Button
        def get_data():
            global Nx, Ny, Lx, Ly, t, dt, steps, nu, x, y, w0
            Nx = int(self.NxEntry.get())
            Ny = int(self.NyEntry.get())
            Lx = check_scalar(self.LxEntry)
            Ly = check_scalar(self.LyEntry)
            t = check_scalar(self.tEntry)
            dt = check_scalar(self.dtEntry)
            steps = int(t / dt)
            nu = check_scalar(self.nuEntry)
            x = np.linspace(-Lx/2, Lx/2, Nx)
            y = np.linspace(-Ly/2, Ly/2, Ny)
            x, y = np.meshgrid(x, y)
            w0 = check_matrix(Ny, Nx, self.w0Entry)
            # Check for unfilled fields
            if None in [Nx, Ny, Lx, Ly, t, dt, steps, nu]:
                return
            if w0 is None:
                return
            # Check for negative values
            pEntries = [self.LxEntry, self.LyEntry, self.tEntry, self.dtEntry]
            for i, positive in enumerate([Lx, Ly, t, dt]):
                if positive <= 0:
                    pEntries[i].delete(0, 'end')
                    self.focus()
            self.quit()

        def check_matrix(i, j, aEntry):
            try:
                a = eval(aEntry.get())
                if not(isinstance(a, np.ndarray)) or np.shape(a)!=(i, j) or a.dtype==np.cfloat:
                    aEntry.delete(0, 'end')
                    self.focus()
                    return 
                return a
            except:
                aEntry.delete(0, 'end')
                self.focus()
                return 

        def check_scalar(aEntry):
            try:
                a = eval(aEntry.get())
                if not(isinstance(a, (int, float))):
                    aEntry.delete(0, 'end')
                    self.focus()
                    return 
                return float(a)
            except:
                aEntry.delete(0, 'end')
                self.focus()
                return 

        self.SimulateButton = ctk.CTkButton(self.frame, text='Simulate', command=get_data)
        self.SimulateButton.grid(row=7, column=2, columnspan=2, padx=20, pady=20, sticky='ew')

def setup():
    ctk.set_appearance_mode('System')  
    ctk.set_default_color_theme('blue') 
    app = App()
    app.mainloop()
    app.destroy()
    return Nx, Ny, Lx, Ly, t, dt, steps, nu, x, y, w0