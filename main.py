import net
import stiffness
import helpers
import gaussNodes
import forceVector
import globalMatrix
import numpy as np
import visualize


xyzLen = [4,2,4]

xyzQuan = [2,1,2]

ZP = [[2,6,-1],[3,6,-1]]

ZU = [[0,5],[1,5]]

npq = xyzQuan[0]*xyzQuan[1]*xyzQuan[2]


AKT = net.createAKT(xyzLen,xyzQuan)
# visualize.plot(AKT)

visualize.plot_AKT_cube(AKT,xyzQuan[0],xyzQuan[1],xyzQuan[2])
NT = net.createNT(AKT)
DFIABG = stiffness.createDFIABG()
jacobiansForElement = [[[[0.0] * 3 for _ in range(3)] for _ in range(npq)] for _ in range(npq)]

for i in range(npq):
    print(i)
    jacobiansForElement[i] = stiffness.createJacobian(NT, AKT, DFIABG, i)


jacobianValuesForElement = [[0.0] * 27 for _ in range(npq)]
for i in range(npq):
    for j in range(27):
        jacobianValuesForElement[i][j] = helpers.getDeterminant(jacobiansForElement[i][j])

DFIXYZ = stiffness.createDFIXYZ(jacobiansForElement, DFIABG, npq)

permutations = gaussNodes.getTwentySevenGaussianNodes()
permutationIndex = {}
for i in range(27):
    permutationIndex[tuple(permutations[i])] = i

allAe11 = stiffness.createAe11(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe22 = stiffness.createAe22(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe33 = stiffness.createAe33(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe12 = stiffness.createAe12(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe13 = stiffness.createAe13(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe23 = stiffness.createAe23(DFIXYZ, jacobianValuesForElement, npq, permutationIndex)

allMGE = stiffness.createMGE(allAe11, allAe22, allAe33, allAe12, allAe13, allAe23, npq)

DPSIET = forceVector.createDPSIET()

permutations2d = gaussNodes.getNineGaussianNodes()
permutationIndex2d = {}
for i in range(9):
    permutationIndex2d[tuple(permutations2d[i])] = i

allFE = [[0.0] * 60 for _ in range(len(ZP))]
for i in range(len(ZP)):
    allFE[i] = forceVector.createFE(DPSIET, AKT, NT, ZP[i], permutationIndex2d)

nqp = len(AKT[0])

# Making all MGE full matrices (tmp to be usable with current gaussian elimination)
for q in range(len(allMGE)):
    for i in range(len(allMGE[0])):
        for j in range(i, len(allMGE[0][0])):
            allMGE[q][j][i] = allMGE[q][i][j]

globalMG = [[0] * (nqp * 3) for _ in range(nqp * 3)]
globalF = [0] * (nqp * 3)

globalMatrix.initGlobals(globalMG, globalF, allMGE, allFE, ZP, NT)
globalMatrix.fixateSides(globalMG, ZU, NT)


coef = 30
print("nice")
print(len(globalMG),len(globalF))

globalMG = np.array(globalMG, dtype=np.float64)
globalF = np.array(globalF, dtype=np.float64)

finalRes = np.linalg.solve(globalMG, globalF)

for i in range(len(AKT[0])):
    AKT[0][i] += coef * finalRes[i * 3]
    AKT[1][i] += coef * finalRes[i * 3 + 1]
    AKT[2][i] += coef * finalRes[i * 3 + 2]
visualize.plot_AKT_cube(AKT)
print("success")


