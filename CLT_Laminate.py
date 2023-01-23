import numpy as np

from CLT_Base import Node
from CLT_Ply import Ply


class Laminate(Node):
    stack = None
    material = None
    thickness = None
    symmetric = None
    repeatLeft = None
    repeatRight = None
    numberOfPlies = None

    Plies = None
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    ABD = np.zeros((6, 6))
    abd = np.zeros((6, 6))
    a = np.zeros((3, 3))
    b = np.zeros((3, 3))
    d = np.zeros((3, 3))

    LengthDisplayOutput = 70
    verbose: bool

    def __init__(self, stack, material, thickness, symmetric, repeatLeft, repeatRight, verbose=False):
        super().__init__(locals())

        self.lam_strain = None

        if self.symmetric:
            self.NumberOfPlies = (len(self.stack) * self.repeatLeft) * 2 * self.repeatRight
        else:
            self.NumberOfPlies = len(self.stack) * self.repeatLeft * self.repeatRight

        if len(self.material) == 1:
            self.material = np.repeat(self.material, len(self.stack))

        if len(self.thickness) == 1:
            self.thickness = np.ones((1, len(self.stack))) * self.thickness

        self.parseStackingSequence()
        self.parsePly()
        self.calculateABD()

    def parseStackingSequence(self):
        stack = self.stack
        repeat_left = self.repeatLeft
        repeat_right = self.repeatRight
        thickness = self.thickness
        material = self.material

        if self.symmetric:
            explicit_stack = np.repeat([np.repeat(stack, repeat_left), np.flip(np.repeat(stack, repeat_left))], repeat_right)
            explicit_thickness = np.repeat(
                [np.repeat(thickness, repeat_left), np.flip(np.repeat(thickness, repeat_left))], repeat_right)
            explicit_material = np.repeat(
                [np.repeat(material, repeat_left), np.flip(np.repeat(material, repeat_left))], repeat_right)
        else:
            # If not symmetric, only repeat right is considered!!
            explicit_stack = np.repeat(stack, repeat_right)
            explicit_thickness = np.repeat(thickness, repeat_right)
            explicit_material = np.repeat(material, repeat_right)

        self.stack = explicit_stack
        self.thickness = explicit_thickness
        self.material = explicit_material

    def parsePly(self):
        plies = []
        for i in range(len(self.stack)):
            plies.append(Ply(angle=self.stack[i],
                             material=self.material[i],
                             thickness=self.thickness[i]))

        self.Plies = plies

    def calculateABD(self):
        # Laminate Local Reference System
        h = sum(self.thickness)
        z = np.zeros((len(self.stack)+1, 1))
        z[0] = -h / 2

        for i in range(1, len(z)):
            z[i] = z[i - 1] + self.thickness[i - 1]

        for i in range(self.NumberOfPlies):
            self.Plies[i].z_MidPlane = z[i] + self.thickness[i] / 2

        for i in range(len(z)-1):
            Qbar = self.Plies[i].Qbar
            self.A = self.A + Qbar * (z[i + 1] - z[i])
            self.B = self.B + (1 / 2) * Qbar * (z[i + 1] ** 2 - z[i] ** 2)
            self.D = self.D + (1 / 3) * Qbar * (z[i + 1] ** 3 - z[i] ** 3)

        for arr in [self.A, self.B, self.D]:
            arr[np.abs(arr) < 1e-8] = 0

        self.ABD = np.vstack((np.hstack((self.A, self.B)),
                              np.hstack((self.B, self.D))))

        self.abd = np.linalg.inv(self.ABD)
        self.a = self.abd[:3, :3]
        self.b = self.abd[:3, 3:]
        self.d = self.abd[3:, 3:]

    def calculateStrain(self, Nx, Ny, Nxy, Mx, My, Mxy):
        F = np.array([Nx, Ny, Nxy, Mx, My, Mxy])
        abd = self.abd.copy()
        abd[np.isnan(abd)] = 0

        self.lam_strain = abd @ np.transpose(F)
        self.lam_strain[np.abs(self.lam_strain) < 1e-14] = 0

        for ply in self.Plies:
            ply.exy = self.lam_strain[:3]
            ply.kxy = self.lam_strain[3:]

    def plyStrains(self):
        h = sum(self.thickness)
        z = np.zeros((len(self.stack) + 1, 1))
        z[0] = -h / 2

        for i in range(1, len(z)):
            z[i] = z[i - 1] + self.thickness[i - 1]
        z = z.flatten()

        eps = self.lam_strain[:3]
        kap = self.lam_strain[3:]
        strains = []
        for zz in z:
            strains.append(eps + zz*kap)
        return strains

    def ProgressiveFailureAnalysis(self, Nx, Ny, Nxy, Mx, My, Mxy, verbose):

        if self.verbose:
            print('Progressive Failure Analysis')
        
        abd = self.abd.copy()
        abd[np.isnan(abd)] = 0
        
        MCD = np.max(np.abs([Nx, Ny, Nxy, Mx, My, Mxy]))
        Nxold = Nx
        F = np.array([Nx, Ny, Nxy, Mx, My, Mxy]) / MCD
        if np.all(F == 0) or np.all(np.isnan(F)):
            print('Initial Loads are unusable. Aborting...\n')
            if self.verbose: print('=' * self.LengthDisplayOutput)
            return

        tmax = MCD
        nsteps = 10000
        dt = np.linspace(0, tmax * 5, nsteps)
        t = 0
        counter = 0
        flag_FPF = False

        test = [ply.FlagFail for ply in self.Plies]
        while not np.all(test) and t < tmax * 5:
            counter += 1
            if counter >= nsteps:
                print('The applied load is not enough to finish the PFA.\nIncrease initial value!\n')
                return
            NM = dt[counter] * F
            
            eps = np.dot(abd, np.transpose(NM))

            for i in range(self.NumberOfPlies):
                self.Plies[i].exy = eps[:3]
                self.Plies[i].kxy = eps[3:]

                self.Plies[i].calculateStress()
                TW = self.material[i].calculateTW(self.Plies[i].S12)
                if TW > 1 and not self.Plies[i].FlagFail:
                    self.Plies[i].failedPly([0.4, 0.4, 0.15])
                    if verbose:
                        print(
                            f'Ply {i} (Angle = {self.Plies[i].angle} deg) FAILED at {(NM[0] / Nxold) * 100}% of given loads!')
                    if not flag_FPF:  # NEEDS TO BE IMPLEMENTED FOR MORE DETAILED FPF
                        print(f'Following First Fly Failure, the laminate would fail at {i}')
                        flag_FPF = True

            test = [ply.FlagFail for ply in self.Plies]

        if verbose:
            print('ALL plies failed!')

        if self.verbose:
            print('Progressive Failure Analysis Successfully Terminated!')
            print('=' * self.LengthDisplayOutput)
    
    def plot_curvature(self):
        import matplotlib.pyplot as plt
        
        x = np.linspace(0, 1, 50)
        y = x.copy()
        w = np.zeros((x.shape[0], x.shape[0]))
        X, Y = np.meshgrid(x, y)

        for i, xx in enumerate(x):
            for j, yy in enumerate(y):
                wx = -self.lam_strain[3] * np.square(xx-np.min(x))
                wy = -self.lam_strain[4] * np.square(yy-np.min(y))
                wxy = -self.lam_strain[5] * (xx-np.min(x)) * (yy-np.min(y))
                w[j, i] = wx + wy + wxy

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(X, Y, w, cmap='jet')
        ax.plot_surface(X, Y, np.zeros(w.shape), alpha=0.25, color='k')
        plt.show()
        
    def plot_strains(self, glob=True, idx=0):
        import matplotlib.pyplot as plt
        plt.figure()
        
        x = []
        y = []
        
        for ply in self.Plies:

            z_top = ply.z_MidPlane - ply.thickness/2
            z_bottom = ply.z_MidPlane + ply.thickness/2
            
            if glob:
                eps_top = ply.exy + z_top * ply.kxy
                eps_top[-1] = eps_top[-1] / 2  # engineering strain
                eps_btm = ply.exy + z_bottom * ply.kxy
                eps_btm[-1] = eps_btm[-1] / 2  # engineering strain
                
                top = np.dot(ply.T, eps_top)
                bottom = np.dot(ply.T, eps_btm)
            else:
                top = ply.exy + z_top * ply.kxy
                bottom = ply.exy + z_bottom * ply.kxy
            
            x = [0, bottom[idx], top[idx], 0]
            y = [z_bottom, z_bottom, z_top, z_top]
            plt.fill(x, y, alpha=0.5, label=str(ply.angle))
        
        ax = plt.gca()
        ax.invert_yaxis()
        plt.legend()
        plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        