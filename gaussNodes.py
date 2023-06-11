import math
import helpers

def getTwentySevenGaussianNodes():
    res = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                res.append([helpers.xg[i], helpers.xg[j], helpers.xg[k]])
    return res

def getNineGaussianNodes():
    res = []
    for i in range(3):
        for j in range(3):
            res.append([helpers.xg[i], helpers.xg[j]])
    return res

