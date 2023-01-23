import numpy as np
from CLT_Base import Node
from CLT_Material import Material


class Ply(Node):
    S12 = None
    Sxy = None
    e12 = None
    exy = None
    kxy = None

    Qbar = None
    Invariants = None
    T = None
    z_MidPlane = None

    FlagFail = False
    Qbar_reduced = None

    angle: float
    material: Material

    def __init__(self, angle, material, thickness):
        super().__init__(locals())

        # self.angle = angle
        # self.material = material
        # self.thickness = thickness

        self.T = self.calculateT(self.angle)
        self.qbar()

    def calculateStress(self):
        self.Sxy = self.Qbar @ (self.exy + self.z_MidPlane * self.kxy)

        self.S12 = self.T @ self.Sxy

        self.e12 = np.dot(self.T, self.exy + self.z_MidPlane * self.kxy)  # This might be wrong if using PFA
        self.e12 = self.T @ (self.exy + self.z_MidPlane * self.kxy)  # This might be wrong if using PFA

    def failedPly(self, Reduction):
        self.FlagFail = True
        self.material.E1 *= Reduction[0]
        self.material.E2 *= Reduction[1]
        self.material.G12 *= Reduction[2]

        self.qbar()

    def qbar(self):
        E1 = self.material.E1
        E2 = self.material.E2
        v12 = self.material.v12
        G12 = self.material.G12
        
        R = np.array([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 2]])
        
        T = self.T.copy()

        v21 = v12 * (E2 / E1)
        q11 = E1 / (1 - (v12 * v21))
        q22 = E2 / (1 - (v12 * v21))
        q12 = (v12 * E2) / (1 - (v12 * v21))
        q66 = G12
        Q = np.array([[q11, q12, 0],
                      [q12, q22, 0],
                      [0, 0, q66]])
        
        self.Qbar = np.linalg.inv(T) @ Q @ R @ T @ np.linalg.inv(R)
        self.Invariants = np.array([(1/8) * (3 * q11 + 3 * q22 + 2 * q12 + 4 * q66),
                                    (1/2) * (q11 - q22),
                                    (1/8) * (q11 + q22 - 2 * q12 - 4 * q66),
                                    (1/8) * (q11 + q22 + 6 * q12 - 4 * q66),
                                    (1/8) * (q11 + q22 - 2 * q12 + 4 * q66)])
    
    @staticmethod
    def calculateT(theta):
        x = np.cos(np.deg2rad(theta))
        y = np.sin(np.deg2rad(theta))

        T = np.array([[x ** 2, y ** 2, 2 * x * y],
                      [y ** 2, x ** 2, -2 * x * y],
                      [-x * y, x * y, x ** 2 - y ** 2]])

        return T
