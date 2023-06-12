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

ZP = [[2,4,-10],[3,2,-10]]

ZU = [[0,5]]


e = 100
nu = 0.3

npq = xyzQuan[0]*xyzQuan[1]*xyzQuan[2]


AKT = net.createAKT(xyzLen,xyzQuan)
# visualize.plot(AKT)
# print(AKT)

# visualize.plot_cube(AKT,xyzLen,xyzQuan)
# visualize.plot_lines_between_dots(AKT,xyzLen,xyzQuan)
NT = net.createNT(AKT,xyzQuan)
# visualize.plot_cube(AKT,NT,xyzQuan)
visualize.plot_cube(AKT,NT,xyzQuan)
print(NT)
# visualize.plot_lines_between_dots(NT,AKT)
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

allAe11 = stiffness.createAe11(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe22 = stiffness.createAe22(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe33 = stiffness.createAe33(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe12 = stiffness.createAe12(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe13 = stiffness.createAe13(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)
allAe23 = stiffness.createAe23(e,nu,DFIXYZ, jacobianValuesForElement, npq, permutationIndex)

allMGE = stiffness.createMGE(allAe11, allAe22, allAe33, allAe12, allAe13, allAe23, npq)

DPSITE = forceVector.createDPSITE()

permutations2d = gaussNodes.getNineGaussianNodes()
permutationIndex2d = {}
for i in range(9):
    permutationIndex2d[tuple(permutations2d[i])] = i

allFE = [[0.0] * 60 for _ in range(len(ZP))]
for i in range(len(ZP)):
    allFE[i] = forceVector.createFE(DPSITE, AKT, NT, ZP[i], permutationIndex2d)

nqp = len(AKT[0])

# fill matrix(mirror elements)
for q in range(len(allMGE)):
    for i in range(len(allMGE[0])):
        for j in range(i, len(allMGE[0][0])):
            allMGE[q][j][i] = allMGE[q][i][j]

globalMG = [[0] * (nqp * 3) for _ in range(nqp * 3)]
globalF = [0] * (nqp * 3)

globalMatrix.initGlobals(globalMG, globalF, allMGE, allFE, ZP, NT)
globalMatrix.fixateSides(globalMG, ZU, NT)


print("nice")
print(len(globalMG),len(globalF))

globalMG = np.array(globalMG, dtype=np.float64)
globalF = np.array(globalF, dtype=np.float64)

U = np.linalg.solve(globalMG, globalF)

for i in range(len(AKT[0])):
    AKT[0][i] += U[i * 3]
    AKT[1][i] +=  U[i * 3 + 1]
    AKT[2][i] += U[i * 3 + 2]
visualize.plot_cube(AKT,NT,xyzQuan)
# print(NT)
# visualize.plot_AKT_cube(AKT,xyzQuan[0],xyzQuan[1],xyzQuan[2])
print("success")


