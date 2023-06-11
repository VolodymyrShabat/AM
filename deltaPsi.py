import helpers

def deltaPsiEta4(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return 0.25 * (1 + tau * taui) * (tau * taui * etai + 2 * eta * pow(etai, 2))

def deltaPsiTau4(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return 0.25 * (1 + eta * etai) * (eta * etai * taui + 2 * tau * pow(taui, 2))

def deltaPsiEta57(eta, tau, nodeLocalNumber):
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return -eta * (1 + tau * taui)

def deltaPsiTau57(eta, tau, nodeLocalNumber):
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return 0.5 * taui * (1 - pow(eta, 2))

def deltaPsiEta68(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    return 0.5 * etai * (1 - pow(tau, 2))

def deltaPsiTau68(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    return -tau * (1 + eta * etai)