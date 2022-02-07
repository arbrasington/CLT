import numpy as np


class Material:
    # stress allowables
    s1t = None
    s1c = None
    s2t = None
    s2c = None
    t12 = None

    # Tsai Wu
    F1 = None
    F11 = None
    F2 = None
    F22 = None
    F12 = None
    F66 = None

    def __init__(self,
                 name,
                 e1,
                 e2,
                 g12,
                 v12,
                 a11=None,
                 a22=None,
                 S=None,
                 Q=None,
                 verbose=False):
        self.E1 = e1
        self.E2 = e2
        self.G12 = g12
        self.v12 = v12
        self.a11 = a11
        self.a22 = a22
        self.S = S
        self.Q = Q
        self.verb = verbose
        self.name = name

    def calculateS(self):
        self.S = np.array([[1 / self.E1, -self.v12 / self.E1, 0],
                           [-self.v12 / self.E1, 1 / self.E2, 0],
                           [0, 0, 1 / self.G12]])

    def calculateQ(self):
        if self.S is None: self.calculateS()

        sub = self.S[0, 0] * self.S[1, 1] - self.S[0, 1] ** 2
        Q11 = self.S[2, 2] / sub
        Q22 = self.S[1, 1] / sub
        Q12 = self.S[1, 2] / sub
        Q66 = 1 / self.S[2, 2]

        self.Q = np.array([[Q11, Q12, 0],
                           [Q12, Q22, 0],
                           [0, 0, Q66]])

    def addFailureTW(self, s1t, s1c, s2t, s2c, t12, F12='auto'):
        self.s1t = s1t
        self.s1c = s1c
        self.s2t = s2t
        self.s2c = s2c
        self.t12 = t12
        self.F12 = F12

        self.F1 = 1 / self.s1t - 1 / self.s1c
        self.F11 = 1 / (self.s1t * self.s1c)
        self.F2 = 1 / self.s2t - 1 / self.s2c
        self.F22 = 1 / (self.s2t * self.s2c)
        self.F66 = 1 / self.t12 ** 2

        if self.F12 == 'auto':
            self.F12 = -(1 / 2) * np.sqrt(1 / (self.s1t * self.s1c * self.s2t * self.s2c))

        if self.F11 * self.F22 - self.F12 ** 2 > 0:
            self.F12 = 0

    def calculateTW(self, S12):
        s1 = S12[0]
        s2 = S12[1]
        t12 = S12[2]

        TW = self.F1 * s1 + self.F2 * s2 + self.F11 * s1 ** 2 + self.F22 * s2 ** 2 + \
             self.F66 * (t12 ** 2) + 2 * self.F12 * s1 * s2

        if TW < 0:
            TW = 0
        return TW
