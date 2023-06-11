import deltaPhi
import gaussNodes
import numpy as np
import helpers


def createDFIABG():
    permutations = gaussNodes.getTwentySevenGaussianNodes()
    DFIABG = [[[0.0] * 20 for _ in range(3)] for _ in range(27)]
    for i in range(27):
        for j in range(8):
            DFIABG[i][0][j] = deltaPhi.deltaPhiAlpha_8(permutations[i][0], permutations[i][1], permutations[i][2], j)
        
        for j in range(8, 20):
            DFIABG[i][0][j] = deltaPhi.deltaPhiAlpha_12(permutations[i][0], permutations[i][1], permutations[i][2], j)
        
        for j in range(8):
            DFIABG[i][1][j] = deltaPhi.deltaPhiBeta_8(permutations[i][0], permutations[i][1], permutations[i][2], j)
        
        for j in range(8, 20):
            DFIABG[i][1][j] = deltaPhi.deltaPhiBeta_12(permutations[i][0], permutations[i][1], permutations[i][2], j)
        
        for j in range(8):
            DFIABG[i][2][j] = deltaPhi.deltaPhiGamma_8(permutations[i][0], permutations[i][1], permutations[i][2], j)
        
        for j in range(8, 20):
            DFIABG[i][2][j] = deltaPhi.deltaPhiGamma_12(permutations[i][0], permutations[i][1], permutations[i][2], j)
    
    return DFIABG


def createJacobian(NT, AKT, DFIABG, element):
    jacobian_matrices_27 = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(27)]
    for node in range(len(DFIABG)):
        for der in range(3):
            for coord in range(3):
                for i in range(20):
                    jacobian_matrices_27[node][der][coord] += AKT[coord][NT[i][element]] * DFIABG[node][der][i]
    return jacobian_matrices_27

def createDFIXYZ(jacobiansForElement, DFIABG, finiteElementsQuantity):
    # q x 27 x 20 x 3
    DFIXYZ = [[[[0.0] * 3 for _ in range(20)] for _ in range(27)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(27):
            for j in range(20):
                DFIXYZ[q][i][j] = np.linalg.solve(jacobiansForElement[q][i], [DFIABG[i][0][j], DFIABG[i][1][j], DFIABG[i][2][j]])
                # print(jacobiansForElement[q][i], " vec b: ", DFIABG[i][0][j], DFIABG[i][1][j], DFIABG[i][2][j], " res: ", DFIXYZ[q][i][j][0], DFIXYZ[q][i][j][1], DFIXYZ[q][i][j][2], '\n')
    return DFIXYZ

def createAe11(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae11 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae11[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] +
                                mu * (DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1] + DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2])
                            ) * jacobianValuesForElement[q][pIndex]
    return ae11


def createAe22(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae22 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae22[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1] +
                                mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] + DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2])
                            ) * jacobianValuesForElement[q][pIndex]
    return ae22


def createAe33(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae33 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae33[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2] +
                                mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] + DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1])
                            ) * jacobianValuesForElement[q][pIndex]
    return ae33


def createAe12(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae12 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae12[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * nu * DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][1] +
                                mu * DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][0]
                            ) * jacobianValuesForElement[q][pIndex]
    return ae12


def createAe13(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae13 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae13[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * nu * DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][2] +
                                mu * DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][0]
                            ) * jacobianValuesForElement[q][pIndex]
    return ae13


def createAe23(DFIXYZ, jacobianValuesForElement, finiteElementsQuantity, permutationIndex):
    e = 100
    nu = 0.3
    lambda_val = e / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    ae23 = [[[0] * 20 for _ in range(20)] for _ in range(finiteElementsQuantity)]
    for q in range(finiteElementsQuantity):
        for i in range(20):
            for j in range(20):
                for m in range(3):
                    for n in range(3):
                        for k in range(3):
                            pIndex = permutationIndex[(helpers.xg[m], helpers.xg[n], helpers.xg[k])]
                            ae23[q][i][j] += helpers.c[m] * helpers.c[n] * helpers.c[k] * (
                                lambda_val * nu * DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][2] +
                                mu * DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][1]
                            ) * jacobianValuesForElement[q][pIndex]
    return ae23

# def createMGE(allAe11, allAe22, allAe33, allAe12, allAe13, allAe23, finiteElementsQuantity):
#     MGE = [[[0] * 60 for _ in range(60)] for _ in range(finiteElementsQuantity)]
#     for q in range(finiteElementsQuantity):
#         for i in range(20):
#             for j in range(20):
#                 MGE[q][i][j] = allAe11[q][i][j]
#                 MGE[q][i + 20][j + 20] = allAe22[q][i][j]
#                 MGE[q][i + 40][j + 40] = allAe33[q][i][j]
#                 MGE[q][i][j + 20] = allAe12[q][i][j]
#                 MGE[q][i][j + 40] = allAe13[q][i][j]
#                 MGE[q][i + 20][j + 40] = allAe23[q][i][j]
#     return MGE


def createMGE(allAe11, allAe22, allAe33, allAe12, allAe13, allAe23, npq):
    MGE = np.zeros((npq, 60, 60))
    MGE[:, :20, :20] = allAe11
    MGE[:, 20:40, 20:40] = allAe22
    MGE[:, 40:, 40:] = allAe33
    MGE[:, :20, 20:40] = allAe12
    MGE[:, :20, 40:] = allAe13
    MGE[:, 20:40, 40:] = allAe23
    return MGE