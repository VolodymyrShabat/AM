import helpers

#F
#x1 y1 z1 x2 y2 z2 ...

#FE
#x1 x2 x3 x4 ... y1 y2 y3 y4 ... z1 z2 z3 z4
def initGlobals(MG, vecF, allMGE, allFE, ZP, NT):
    for q in range(len(allMGE)):
        for i in range(len(allMGE[0])):
            for j in range(len(allMGE[0])):
                i1 = NT[i % 20][q] * 3 + i // 20 #[i % 20] - move xn,yn,zn  i // 20 - move between x,y,z
                i2 = NT[j % 20][q] * 3 + j // 20 # same but for column
                MG[i1][i2] += allMGE[q][i][j]

    for q in range(len(ZP)):
        print(len(allFE[0]))
        for i in range(len(allFE[0])):
            i1 = NT[i % 20][ZP[q][0]] * 3 +  i // 20
            vecF[i1] += allFE[q][i]


#MG
#x1 x2 x3
#x2
#x3

#MG
#x1 y1 z1
#y1
#z1
def fixateSides(globalMatrix, ZU, NT):
    for i in range(len(ZU)):
        for j in range(8):
            index = NT[helpers.sideToPoint[ZU[i][1]][j]][ZU[i][0]] * 3
            globalMatrix[index][index] = 10 ** 30
            globalMatrix[index + 1][index + 1] = 10 ** 30
            globalMatrix[index + 2][index + 2] = 10 ** 30

