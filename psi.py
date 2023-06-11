import helpers

def psi4(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return 0.25 * (1 + eta * etai) * (1 + tau * taui) * (eta * etai + tau * taui - 1)


def psi57(eta, tau, nodeLocalNumber):
    taui = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][1]
    return 0.5 * (1 - eta**2) * (1 + tau * taui)


def psi68(eta, tau, nodeLocalNumber):
    etai = helpers.nodeNumberToLocalCoords2d[nodeLocalNumber][0]
    return 0.5 * (1 - tau**2) * (1 + eta * etai)


def psii(eta, tau, nodeLocalNumber, i):
    if i in [0, 1, 2, 3]:
        return psi4(eta, tau, nodeLocalNumber)
    elif i in [4, 6]:
        return psi57(eta, tau, nodeLocalNumber)
    elif i in [5, 7]:
        return psi68(eta, tau, nodeLocalNumber)