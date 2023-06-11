def createAKT(xyzLength, xyzQuantity):
    xyzDistance = [xyzLength[0] / xyzQuantity[0], xyzLength[1] / xyzQuantity[1], xyzLength[2] / xyzQuantity[2]]
    result = [[], [], []]
    layersQuantity = 1 + 2 * xyzQuantity[2]
    xyzDistance[2] = xyzLength[2] / (layersQuantity - 1)
    xyzDistance[1] = xyzLength[1] / ((1 + 2 * xyzQuantity[1]) - 1)
    xyzDistance[0] = xyzLength[0] / ((1 + 2 * xyzQuantity[0]) - 1)

    for layer in range(layersQuantity):
        if layer % 2 == 0:
            for row in range(1 + 2 * xyzQuantity[1]):
                if row % 2 == 0:
                    for i in range(1 + 2 * xyzQuantity[0]):
                        result[0].append(i * xyzDistance[0])
                        result[1].append(row * xyzDistance[1])
                        result[2].append(layer * xyzDistance[2])
                else:
                    for i in range(1 + xyzQuantity[0]):
                        result[0].append(i * xyzDistance[0] * 2)
                        result[1].append(row * xyzDistance[1])
                        result[2].append(layer * xyzDistance[2])
        else:
            for row in range(1 + xyzQuantity[1]):
                for i in range(1 + xyzQuantity[0]):
                    result[0].append(i * xyzDistance[0] * 2)
                    result[1].append(row * xyzDistance[1] * 2)
                    result[2].append(layer * xyzDistance[2])
    
    return result

def createNT(AKT, xyzQuantity):
    totalQuantity = xyzQuantity[0] * xyzQuantity[1] * xyzQuantity[2]
    longRowPointsCount = 1 + 2 * xyzQuantity[0]
    shortRowPointsCount = 1 + xyzQuantity[0]
    bigLayerPointsCount = ((longRowPointsCount + shortRowPointsCount) * xyzQuantity[1]) + longRowPointsCount
    smallLayerPointsCount = shortRowPointsCount * (xyzQuantity[1] + 1)
    bothLayersPointsCount = bigLayerPointsCount + smallLayerPointsCount
    bothRowsPointsCount = longRowPointsCount + shortRowPointsCount
    NT = [[0] * totalQuantity for _ in range(20)]
    cubesInLayer = xyzQuantity[0] * xyzQuantity[1]

    for z in range(xyzQuantity[2]):
        for y in range(xyzQuantity[1]):
            for x in range(xyzQuantity[0]):
                cubeNumber = cubesInLayer * z + xyzQuantity[0] * y + x
                NT[0][cubeNumber] = z * bothLayersPointsCount + y * bothRowsPointsCount + x * 2
                NT[1][cubeNumber] = NT[0][cubeNumber] + 2
                NT[2][cubeNumber] = NT[1][cubeNumber] + bothRowsPointsCount
                NT[3][cubeNumber] = NT[2][cubeNumber] - 2

                for i in range(4, 8):
                    NT[i][cubeNumber] = NT[i - 4][cubeNumber] + bothLayersPointsCount

                NT[8][cubeNumber] = NT[0][cubeNumber] + 1
                NT[9][cubeNumber] = z * bothLayersPointsCount + y * bothRowsPointsCount + longRowPointsCount + x + 1
                NT[10][cubeNumber] = NT[8][cubeNumber] + bothRowsPointsCount
                NT[11][cubeNumber] = NT[9][cubeNumber] - 1

                NT[12][cubeNumber] = z * bothLayersPointsCount + bigLayerPointsCount + y * shortRowPointsCount + x
                NT[13][cubeNumber] = NT[12][cubeNumber] + 1
                NT[14][cubeNumber] = NT[13][cubeNumber] + shortRowPointsCount
                NT[15][cubeNumber] = NT[14][cubeNumber] - 1

                for i in range(16, 20):
                    NT[i][cubeNumber] = NT[i - 8][cubeNumber] + bothLayersPointsCount
    
    return NT