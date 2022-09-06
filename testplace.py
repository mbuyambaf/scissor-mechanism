import numpy as np
import csv

# check that member is in equalibrium


def check(*args):
    checklist = list()
    for i in args:
        checklist.append(i)

    total = round(sum(checklist), 10)

    if total != 0:
        print(f'This member is unstable! Resultant = {total}')
    else:
        print(f'This member is stable! Resultant = {total}')


# P
# physical constants
# material = 1  # mass(kg) per length
#g = 9.81
g = 9.81
g1 = 9.81
rho = 8000  # density if stainless steel AISI 316 Stainless Steel (kg/m3)

# geometric propteries
# AB
l_AB = 1.28
vol_AB = 0.0022825499
# rails (x2), shafts (x3)
# additionPlatformMass = ((0.000288073749*2) + (0.00008123*3))*rho
additionPlatformMass = (0)*rho

m_AB = vol_AB*rho + additionPlatformMass
w_AB = m_AB*g


# CD
l_CD = 1.000
vol_CD = 0.000394883597
m_CD = vol_CD*rho*2
w_CD = m_CD*g

# EB
l_EB = 1.000
#m_EB = material*l_EB
vol_EB = 0.000394883597
m_EB = vol_EB*rho*2
w_EB = m_EB*g


##############################################################################
# inputs
positions = 2

#startAngle = 2.48*np.pi/180
startAngle = 13.907*np.pi/180
endAngle = 46.58*np.pi/180
thetaList = np.linspace(endAngle, startAngle, num=positions)

gidMass = 0.350  # grams
fullBin = 252
FOS = 1
PList = np.linspace(gidMass*g1, gidMass*fullBin*g1*FOS, num=positions)

L1 = 0.300101  # connecting link length
# L1 = 0.400101  # connecting link length
# L2 = 0.148592  # fixing point from D
L2 = 0.15  # fixing point from D
L3 = 0.15  # veritical displacement


print("Results! \n")
for theta, P in zip(thetaList, PList):
    # secondary variables
    L4 = (l_CD/2)-L2
    alpha = np.arcsin((L3 + L2*np.sin(theta))/L1)

    # calculation
    # Link B
    Rbx = 0
    Rby = (0.5*l_AB*(P+w_AB)-(l_CD*np.cos(theta)*(P+w_AB)))/(l_CD*np.cos(theta))
    Rcx = 0
    Rcy = w_AB + P - Rby
    check(Rby, -P, -w_AB, Rcy, Rbx, Rcx)

    # Link EB
    Rfx = 0
    Rfy = -Rby - w_EB
    Rex = 0
    Rey = Rfy + w_EB + Rby
    check(-Rbx, -Rfx, Rex, -w_EB, Rey, -Rfy, -Rby)

    # Link CD
    F_hp = (0.5*l_CD*(2*Rcy - Rfy + w_CD))/(-0.5*l_CD*np.sin(alpha) +
                                            (L4-0.5*l_CD)*np.tan(theta)*np.cos(alpha) + L4*np.sin(alpha))
    Rdx = - F_hp * np.cos(alpha)
    Rdy = Rcy - Rfy + w_CD - F_hp*np.sin(alpha)
    check(Rdx, F_hp*np.cos(alpha), Rfx, -Rcy+Rfy-w_CD, Rdy+F_hp*np.sin(alpha))

    #results (outputs)
    requiredForce = round(F_hp*np.cos(alpha), 2)
    angularPosition = round(theta*180/np.pi, 2)
    massOfGIDs = round(P/g1, 2)
    numberGID = round(massOfGIDs/gidMass, 0)

    print(
        f'Required pump force: {requiredForce} at angle {angularPosition} degrees with mass {massOfGIDs} kg')
    requiredForce2 = round(F_hp*np.sin(alpha), 2)
    # print(requiredForce2)
    # print(alpha*180/np.pi)
    if abs(requiredForce) <= 6000:
        print('This is the ideal position!!!!!!!!!!!!!!!!!!')
