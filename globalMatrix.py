import helpers

def initGlobals(globalMatrix, vecF, allMGE, allFE, ZP, NT):
    for q in range(len(allMGE)):
        for i in range(len(allMGE[0])):
            for j in range(len(allMGE[0])):
                compx = i // 20
                compy = j // 20
                i1 = NT[i % 20][q] * 3 + compx
                i2 = NT[j % 20][q] * 3 + compy
                globalMatrix[i1][i2] += allMGE[q][i][j]

    for q in range(len(ZP)):
        for i in range(len(allFE[0])):
            comp = i // 20
            i1 = NT[i % 20][ZP[q][0]] * 3 + comp
            vecF[i1] += allFE[q][i]

def fixateSides(globalMatrix, ZU, NT):
    for i in range(len(ZU)):
        for j in range(8):
            index = NT[helpers.sideToPoint[ZU[i][1]][j]][ZU[i][0]] * 3

            globalMatrix[index][index] = 10 ** 30
            globalMatrix[index + 1][index + 1] = 10 ** 30
            globalMatrix[index + 2][index + 2] = 10 ** 30