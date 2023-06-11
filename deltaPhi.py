import helpers
def deltaPhiAlpha_8(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.125 * (1 + beta * bi) * (1 + gamma * gi) * (ai * (2 * alpha * ai + beta * bi + gamma * gi - 1))

def deltaPhiAlpha_12(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.25 * (1 + beta * bi) * (1 + gamma * gi) * (ai - 2 * alpha * bi**2 * gi**2 - 3 * alpha**2 * ai * bi**2 * gi**2 - ai * (beta * ai * gi)**2 - ai * (gamma * ai * bi)**2)

def deltaPhiBeta_8(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.125 * (1 + alpha * ai) * (1 + gamma * gi) * (bi * (alpha * ai + 2 * beta * bi + gamma * gi - 1))

def deltaPhiBeta_12(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.25 * (1 + alpha * ai) * (1 + gamma * gi) * (bi - bi * (alpha * bi * gi)**2 - 2 * beta * ai**2 * gi**2 - 3 * beta**2 * ai**2 * bi * gi**2 - bi * (gamma * ai * bi)**2)

def deltaPhiGamma_8(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.125 * (1 + alpha * ai) * (1 + beta * bi) * (gi * (alpha * ai + beta * bi + 2 * gamma * gi - 1))

def deltaPhiGamma_12(alpha, beta, gamma, nodeLocalNumber):
    coords = helpers.nodeNumberToLocalCoords[nodeLocalNumber]
    ai, bi, gi = coords[0], coords[1], coords[2]
    return 0.25 * (1 + alpha * ai) * (1 + beta * bi) * (gi * (1 - (alpha * bi * gi)**2 - (beta * ai * gi)**2) - 2 * gamma * ai**2 * bi**2 - 3 * gamma**2 * gi * ai**2 * bi**2)
