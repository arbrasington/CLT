import numpy as np


class Ply:
    S12 = None
    Sxy = None
    e12 = None
    exy = None
    kxy = None

    Qbar = None
    T = None
    z_MidPlane = None

    FlagFail = False
    Qbar_reduced = None

    def __init__(self,
                 angle,
                 material,
                 thickness):
        self.angle = angle
        self.material = material
        self.thickness = thickness

        self.Qbar = self.qbar(self.material.E1,
                              self.material.E2,
                              self.material.v12,
                              self.material.G12,
                              self.angle)
        self.T = self.calculateT(self.angle)

    def calculateStress(self):
        self.Sxy = np.dot(self.Qbar, self.exy + self.z_MidPlane * self.kxy)

        self.S12 = np.dot(self.T, self.Sxy)
        self.e12 = np.dot(self.T, self.exy + self.z_MidPlane * self.kxy)  # This might be wrong if using PFA
        # self.e12 = np.dot(self.material.S, self.S12)  # This might be wrong if using PFA

    def failedPly(self, Reduction):
        self.FlagFail = True
        self.Qbar = self.qbar(self.material.E1 * Reduction[0],
                              self.material.E2 * Reduction[1],
                              self.material.v12,
                              self.material.G12 * Reduction[2],
                              self.angle)

    @staticmethod
    def qbar(E1, E2, v12, G12, theta):
        v21 = v12 * (E2 / E1)
        q11 = E1 / (1 - (v12 * v21))
        q22 = E2 / (1 - (v12 * v21))
        q12 = (v12 * E2) / (1 - (v12 * v21))
        q66 = G12

        c = np.cos(np.deg2rad(theta))
        s = np.sin(np.deg2rad(theta))

        qb11 = (q11 * c ** 4) + (2 * (q12 + 2 * q66) * s ** 2 * c ** 2) + (q22 * s ** 4)
        qb12 = ((q11 + q22 - 4 * q66) * s ** 2 * c ** 2) + (q12 * (s ** 4 + c ** 4))
        qb22 = (q11 * s ** 4) + (2 * (q12 + 2 * q66) * s ** 2 * c ** 2) + (q22 * c ** 4)
        qb16 = ((q11 - q12 - 2 * q66) * s * c ** 3) + ((q12 - q22 + 2 * q66) * s ** 3 * c)
        qb26 = ((q11 - q12 - 2 * q66) * s ** 3 * c) + ((q12 - q22 + 2 * q66) * s * c ** 3)
        qb66 = ((q11 + q22 - 2 * q12 - 2 * q66) * s ** 2 * c ** 2) + (q66 * (s ** 4 + c ** 4))

        qb = np.array([[qb11, qb12, qb16],
                       [qb12, qb22, qb26],
                       [qb16, qb26, qb66]])

        return qb

    @staticmethod
    def calculateT(theta):
        x = np.cos(np.deg2rad(theta))
        y = np.sin(np.deg2rad(theta))

        T = np.array([[x ** 2, y ** 2, 2 * x * y],
                      [y ** 2, x ** 2, -2 * x * y],
                      [-x * y, x * y, x ** 2 - y ** 2]])

        return T
