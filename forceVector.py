import numpy as np
import gaussNodes
import deltaPsi
import helpers
import psi

def createDPSIET():
    nineNodes = gaussNodes.getNineGaussianNodes()
    DPSIET = np.zeros((9, 2, 8))
    for i in range(9):
        for j in range(4):
            DPSIET[i][0][j] = deltaPsi.deltaPsiEta4(nineNodes[i][0], nineNodes[i][1], j)
        DPSIET[i][0][4] = deltaPsi.deltaPsiEta57(nineNodes[i][0], nineNodes[i][1], 4)
        DPSIET[i][0][6] = deltaPsi.deltaPsiEta57(nineNodes[i][0], nineNodes[i][1], 6)
        DPSIET[i][0][5] = deltaPsi.deltaPsiEta68(nineNodes[i][0], nineNodes[i][1], 5)
        DPSIET[i][0][7] = deltaPsi.deltaPsiEta68(nineNodes[i][0], nineNodes[i][1], 7)
        for j in range(4):
            DPSIET[i][1][j] = deltaPsi.deltaPsiTau4(nineNodes[i][0], nineNodes[i][1], j)
        DPSIET[i][1][4] = deltaPsi.deltaPsiTau57(nineNodes[i][0], nineNodes[i][1], 4)
        DPSIET[i][1][6] = deltaPsi.deltaPsiTau57(nineNodes[i][0], nineNodes[i][1], 6)
        DPSIET[i][1][5] = deltaPsi.deltaPsiTau68(nineNodes[i][0], nineNodes[i][1], 5)
        DPSIET[i][1][7] = deltaPsi.deltaPsiTau68(nineNodes[i][0], nineNodes[i][1], 7)
    return DPSIET

def createDPSIXYZ(DPSIET, AKT, NT, finiteElementIndex, side):
    DPSIXYZ = np.zeros((9, 2, 3))
    for i in range(9):
        for j in range(2):
            for k in range(3):
                for s in range(8):
                    DPSIXYZ[i][j][k] += AKT[k][NT[helpers.sideToPoint[side][s]][finiteElementIndex]] * DPSIET[i][j][s]
    return DPSIXYZ

def createFe1(DPSIET, AKT, NT, ZPi, permutationIndex):
    Pn = ZPi[2]
    fe1 = [0] * 8
    DPSIXYZ = createDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1])
    for i in range(8):
        for m in range(3):
            for n in range(3):
                pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n])]
                nul = DPSIXYZ[pIndex][0][1] * DPSIXYZ[pIndex][1][2] - DPSIXYZ[pIndex][0][2] * DPSIXYZ[pIndex][1][1]
                fe1[i] += helpers.c[m] * helpers.c[n] * Pn * nul * psi.psii(helpers.xg[m], helpers.xg[n], i, i)
    return fe1


def createFe2(DPSIET, AKT, NT, ZPi, permutationIndex):
    Pn = ZPi[2]
    fe2 = [0] * 8
    DPSIXYZ = createDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1])
    for i in range(8):
        for m in range(3):
            for n in range(3):
                pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n])]
                fe2[i] += helpers.c[m] * helpers.c[n] * Pn * (DPSIXYZ[pIndex][0][2] * DPSIXYZ[pIndex][1][0] - DPSIXYZ[pIndex][0][0] * DPSIXYZ[pIndex][1][2]) * psi.psii(helpers.xg[m], helpers.xg[n], i, i)
    return fe2


def createFe3(DPSIET, AKT, NT, ZPi, permutationIndex):
    Pn = ZPi[2]
    fe3 = [0] * 8
    DPSIXYZ = createDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1])
    for i in range(8):
        for m in range(3):
            for n in range(3):
                pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n])]
                fe3[i] += helpers.c[m] * helpers.c[n] * Pn * (DPSIXYZ[pIndex][0][0] * DPSIXYZ[pIndex][1][1] - DPSIXYZ[pIndex][0][1] * DPSIXYZ[pIndex][1][0]) * psi.psii(helpers.xg[m], helpers.xg[n], i, i)
    return fe3


def createFE(DPSIET, AKT, NT, ZPi, permutationIndex):
    fe1 = createFe1(DPSIET, AKT, NT, ZPi, permutationIndex)
    fe2 = createFe2(DPSIET, AKT, NT, ZPi, permutationIndex)
    fe3 = createFe3(DPSIET, AKT, NT, ZPi, permutationIndex)
    fe = [0] * 60
    for i in range(8):
        fe[helpers.sideToPoint[ZPi[1]][i]] = fe1[i]
        fe[helpers.sideToPoint[ZPi[1]][i] + 20] = fe2[i]
        fe[helpers.sideToPoint[ZPi[1]][i] + 40] = fe3[i]
    return fe