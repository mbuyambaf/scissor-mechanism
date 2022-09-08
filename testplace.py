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


# physical constants
# material = 1  # mass(kg) per length
g = 9.81
g1 = 9.81
rho = 8000  # density if stainless steel AISI 316 Stainless Steel (kg/m3)

# geometric propteries
# AB
l_AB = 1.285
vol_AB = 0.0022825499

# rails (x2), shafts (x3)
additionPlatformMass = ((0.000288073749*2) + (0.00008123*3))*rho
#additionPlatformMass = (0)*rho

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
theta2 = 13.3*np.pi/180
theta1 = 2.63365*np.pi/180

gidMass = 0.350  # grams
fullBin = 252
FOS = 1
P = gidMass*fullBin*g1*FOS

L1 = 0.400101  # connecting link length
L2 = 0.148592  # fixing point from D
L3 = 0.127  # veritical displacement


print("Results! \n")

# secondary variables
L4 = (l_CD/2)-L2
alpha = np.arcsin((L3 + L2*np.sin(theta2))/L1)

a = np.matrix(
    [
        [1, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0],
        [0, -1, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, np.cos(alpha)],
        [0, 0, -1, 0, 0, 1, np.sin(alpha)],
        [0, 0, l_CD*np.cos(theta1), 0, 0, 0, L2*(np.sin(alpha)
                                                 * np.cos(theta2)+np.sin(theta2)*np.cos(alpha))],
    ]
)
b = np.matrix(
    [
        [0],
        [P+w_AB],
        [0.5*l_AB*(P+w_AB)/(l_CD*np.cos(theta1))],
        [w_EB],
        [0],
        [w_CD],
        [-0.5*w_CD*l_CD*np.cos(theta1)]
    ]
)

x = np.linalg.solve(a, b)
print(x)

requiredForcex = round(x[6, 0]*np.cos(alpha), 2)
requiredForcey = round(x[6, 0]*np.sin(alpha), 2)
angularPosition = round(theta2*180/np.pi, 2)
massOfGIDs = round(P/g1, 2)
numberGID = round(massOfGIDs/gidMass, 0)

print(
    f'Required pump force: {requiredForcex} at angle {angularPosition} degrees with mass {massOfGIDs} kg')
print(requiredForcex)
print(requiredForcey)
# print(alpha*180/np.pi)
if abs(requiredForcex) <= 1000*g1/FOS:
    print('This Meets your Saferty Factor')
