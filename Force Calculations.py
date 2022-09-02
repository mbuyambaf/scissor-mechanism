import numpy as np
import csv 

data = list()
header = ['Required Force','Angle position','Number of GIDs','Mass of GIDs']

#check that member is in equalibrium
def check(*args):
    checklist = list()
    for i in args:
        checklist.append(i)

    total = round(sum(checklist),10)

    if total != 0:
        print(f'This member is unstable! Resultant = {total}')
    else: 
        print(f'This member is stable! Resultant = {total}')
    

##############################################################################P
# physical constants
#material = 1  # mass(kg) per length
g = 9.81
rho = 8000 # density if stainless steel AISI 316 Stainless Steel (kg/m3)
#theta = 46.58*np.pi/180

##############################################################################P
# geometric propteries
# AB
l_AB = 1.28
vol_AB = 0.0022825499
additionPlatformMass = ((0.000288073749*2) + (0.00008123*3))*rho # rails (x2), shafts (x3)

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
startAngle = 32.345*np.pi/180
endAngle = 46.58*np.pi/180
thetaList = np.linspace(endAngle, startAngle, num=positions)

gidMass = 0.350 #grams
fullBin = 252
FOS = 2
PList = np.linspace(gidMass*g,gidMass*fullBin*g*FOS,num=positions)

# L1 = 0.58362 # connecting link length
# L2 = 0.12554 # fixing point from D
# L3 = 0.095352 # veritical displacement

L1 = 0.4 # connecting link length
L2 = 0.148 # fixing point from D
L3 = 0.163 # veritical displacement


print("Results! \n")
for theta,P in zip(thetaList, PList):
    # secondary variables
    L4 = (l_CD/2)-L2
    alpha = np.arcsin((L3 + L2*np.sin(theta))/L1)
    
    # calculation
    #Link B
    Rbx = 0
    Rby = (0.5*l_AB*(P+w_AB)-(l_CD*np.cos(theta)*(P+w_AB)))/(l_CD*np.cos(theta))
    Rcx = 0
    Rcy = w_AB + P - Rby
    check(Rby,-P,-w_AB,Rcy,Rbx,Rcx)

    # Link EB
    Rfx = 0
    Rfy = -2*Rby - w_EB
    Rex = 0
    Rey = Rfy + w_EB + Rby 
    check(-Rbx,-Rfx,Rex,-w_EB,Rey,-Rfy,-Rby)

    #Link CD
    F_hp = (0.5*l_CD*(2*Rcy - Rfy + w_CD))/(-0.5*l_CD*np.sin(alpha) + (L4-0.5*l_CD)*np.tan(theta)*np.cos(alpha) + L4*np.sin(alpha))
    Rdx = - F_hp *np.cos(alpha)
    Rdy = Rcy - Rfy + w_CD -F_hp*np.sin(alpha)
    check(Rdx,F_hp*np.cos(alpha),Rfx,-Rcy+Rfy-w_CD,Rdy+F_hp*np.sin(alpha))

    #results (outputs)
    requiredForce = round(F_hp*np.cos(alpha),2)
    angularPosition = round(theta*180/np.pi,2)
    massOfGIDs = round(P/g,2)
    numberGID = round(massOfGIDs/gidMass,0)
    

    print(f'Required pump force: {requiredForce} at angle {angularPosition} degrees with mass {massOfGIDs} kg')
    requiredForce2 = round(F_hp*np.sin(alpha),2)
    print(requiredForce2)
    #print(F_hp)
    #rowData = [requiredForce,angularPosition,numberGID,massOfGIDs]
    #data.append(rowData)

    #newLimit = abs(requiredForce)
    
    #print(f'{counter}-New position for L2 is {round(L2,5)} for L3 is {round(L3,5)}. New pump force: {newLimit}')
    #L2 = L2 + 0.001
    #L3 = L3 + 0.001
    #counter+=1


# if __name__ == "__main__":
#     main()
# with open('Bin data.csv','w', encoding='UTF8', newline='') as f:
#     writer = csv.writer(f)

#     writer.writerow(header)

#     writer.writerows(data)


#print(data)