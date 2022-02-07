import numpy as np

from CLT_Ply import Ply


class Laminate:
    Stack = None
    Material = None
    Thickness = None
    Symmetric = None
    RepeatLeft = None
    RepeatRight = None
    NumberOfPlies = None

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
    Verbose = False

    def __init__(self,
                 stack,
                 material,
                 thickness,
                 symmetric,
                 repeat_left,
                 repeat_right,
                 verbose=False):
        self.Stack = stack
        self.Material = material
        self.Thickness = thickness
        self.Symmetric = symmetric
        self.RepeatLeft = repeat_left
        self.RepeatRight = repeat_right
        self.Verbose = verbose

        self.lam_strain = None

        if self.Symmetric:
            self.NumberOfPlies = (len(self.Stack) * self.RepeatLeft) * 2 * self.RepeatRight
        else:
            self.NumberOfPlies = len(self.Stack) * self.RepeatLeft * self.RepeatRight

        if len(self.Material) == 1:
            self.Material = np.repeat(self.Material, len(self.Stack))

        if len(self.Thickness) == 1:
            self.Thickness = np.ones((1, len(self.Stack))) * self.Thickness

        self.parseStackingSequence()
        self.parsePly()
        self.calculateABD()

    def parseStackingSequence(self):
        stack = self.Stack
        repeat_left = self.RepeatLeft
        repeat_right = self.RepeatRight
        thickness = self.Thickness
        material = self.Material

        if self.Symmetric:
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

        self.Stack = explicit_stack
        self.Thickness = explicit_thickness
        self.Material = explicit_material

    def parsePly(self):
        plies = []
        for i in range(len(self.Stack)):
            plies.append(Ply(angle=self.Stack[i],
                             material=self.Material[i],
                             thickness=self.Thickness[i]))

        self.Plies = plies

    def calculateABD(self):
        # Laminate Local Rerefence System
        h = sum(self.Thickness)
        z = np.zeros((len(self.Stack)+1, 1))
        z[0] = -h / 2

        for i in range(1, len(z)):
            z[i] = z[i - 1] + self.Thickness[i - 1]

        for i in range(self.NumberOfPlies):
            self.Plies[i].z_MidPlane = z[i] + self.Thickness[i] / 2

        for i in range(len(z)-1):
            Qbar = self.Plies[i].Qbar
            self.A = self.A + Qbar * (z[i + 1] - z[i])
            self.B = self.B + (1 / 2) * Qbar * (z[i + 1] ** 2 - z[i] ** 2)
            self.D = self.D + (1 / 3) * Qbar * (z[i + 1] ** 3 - z[i] ** 3)

        self.a = np.linalg.inv(self.A)
        if np.all(abs(self.B) <= 1e-6):
            self.b = np.empty((3, 3))
            self.b.fill(np.nan)
        else:
            self.b = np.linalg.inv(self.B)

        self.d = np.linalg.inv(self.D)

        self.ABD = np.vstack((np.hstack((self.A, self.B)),
                              np.hstack((self.B, self.D))))
        # self.abd = np.vstack((np.hstack((self.a, self.b)),
        #                       np.hstack((self.b, self.d))))
        self.abd = np.linalg.inv(self.ABD)
        self.a = self.abd[:3, :3]
        self.b = self.abd[:3, 3:]
        self.d = self.abd[3:, 3:]

    def calculateStrain(self, Nx, Ny, Nxy, Mx, My, Mxy):
        F = np.array([Nx, Ny, Nxy, Mx, My, Mxy])
        abd = self.abd.copy()
        abd[np.isnan(abd)] = 0
        self.lam_strain = np.dot(abd, np.transpose(F))

    def plyStrains(self):
        h = sum(self.Thickness)
        z = np.zeros((len(self.Stack) + 1, 1))
        z[0] = -h / 2

        for i in range(1, len(z)):
            z[i] = z[i - 1] + self.Thickness[i - 1]
        z = z.flatten()

        eps = self.lam_strain[:3]
        kap = self.lam_strain[3:]
        strains = []
        for zz in z:
            strains.append(eps + zz*kap)
        print(strains)
        # return

    def ProgressiveFailureAnalysis(self, Nx, Ny, Nxy, Mx, My, Mxy, verbose):

        if self.Verbose:
            print('Progressive Failure Analysis')

        MCD = np.max(np.abs([Nx, Ny, Nxy, Mx, My, Mxy]))
        Nxold = Nx
        F = np.array([Nx, Ny, Nxy, Mx, My, Mxy]) / MCD
        if np.all(F == 0) or np.all(np.isnan(F)):
            print('Initial Loads are unusable. Aborting...\n')
            if self.Verbose: print('=' * self.LengthDisplayOutput)
            return

        tmax = MCD
        nsteps = 10000
        dt = np.linspace(0, tmax * 5, nsteps)
        t = 0
        counter = 0
        flag_FPF = False

        test = [ply.FlagFail for ply in self.Plies]
        while not np.all(test) and t < tmax * 5:
            counter = counter + 1
            if counter > nsteps:
                print('The applied load is not enough to finish the PFA.\nIncrease initial value!\n')
                return
            NM = dt[counter] * F
            eps, _, _, _ = np.linalg.lstsq(self.ABD, np.transpose(NM), rcond=None)

            for i in range(self.NumberOfPlies):
                self.Plies[i].exy = eps[:3]  # TODO: might need to change the end +=1
                self.Plies[i].kxy = eps[3:]

                self.Plies[i].calculateStress()
                TW = self.Material[i].calculateTW(self.Plies[i].S12)
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

        if self.Verbose:
            print('Progressive Failure Analysis Successfully Terminated!')
            print('=' * self.LengthDisplayOutput)
