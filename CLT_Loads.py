import numpy as np

from CLT_Laminate import Laminate
from CLT_Material import Material


def thermal_loads(dT, laminate: Laminate):
    NxxT = 0
    NyyT = 0
    NxyT = 0
    MxxT = 0
    MyyT = 0
    MxyT = 0

    h = sum(laminate.Thickness)
    z = np.zeros((len(laminate.Stack) + 1, 1))
    z[0] = -h / 2

    for i in range(1, len(z)):
        z[i] = z[i - 1] + laminate.Thickness[i - 1]
    z = z.flatten()

    for i, ply in enumerate(laminate.Plies):
        mat: Material = ply.material
        axx = mat.a11 * np.cos(np.deg2rad(ply.angle)) ** 2 + \
              mat.a22 * np.sin(np.deg2rad(ply.angle)) ** 2
        ayy = mat.a11 * np.sin(np.deg2rad(ply.angle)) ** 2 + \
              mat.a22 * np.cos(np.deg2rad(ply.angle)) ** 2
        axy = 2 * np.cos(np.deg2rad(ply.angle)) * np.sin(np.deg2rad(ply.angle)) * (mat.a11 - mat.a22)

        nxxt = (ply.Qbar[0, 0] * axx + ply.Qbar[0, 1] * ayy + ply.Qbar[0, 2] * axy) * (z[i + 1] - z[i])
        nyyt = (ply.Qbar[1, 0] * axx + ply.Qbar[1, 1] * ayy + ply.Qbar[1, 2] * axy) * (z[i + 1] - z[i])
        nxyt = (ply.Qbar[2, 0] * axx + ply.Qbar[2, 1] * ayy + ply.Qbar[2, 2] * axy) * (z[i + 1] - z[i])

        mxxt = (ply.Qbar[0, 0] * axx + ply.Qbar[0, 1] * ayy + ply.Qbar[0, 2] * axy) * (z[i + 1]**2 - z[i]**2)
        myyt = (ply.Qbar[1, 0] * axx + ply.Qbar[1, 1] * ayy + ply.Qbar[1, 2] * axy) * (z[i + 1]**2 - z[i]**2)
        mxyt = (ply.Qbar[2, 0] * axx + ply.Qbar[2, 1] * ayy + ply.Qbar[2, 2] * axy) * (z[i + 1]**2 - z[i]**2)

        NxxT += nxxt
        NyyT += nyyt
        NxyT += nxyt
        MxxT += mxxt
        MyyT += myyt
        MxyT += mxyt
    return NxxT*dT, NyyT*dT, NxyT*dT, MxxT*dT/2, MyyT*dT/2, MxyT*dT/2
