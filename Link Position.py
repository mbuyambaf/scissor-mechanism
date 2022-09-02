from ast import arguments
from cmath import nan
import numpy as np
import csv 
import math
from tqdm import tqdm

#check that member is in equalibrium
def check(*args):
    checklist = list()
    for i in args:
        checklist.append(i)

    total = round(sum(checklist),5)

    if total != 0:
        print(f'This member is unstable! Resultant = {total}')
    else: 
        print(f'This member is stable! Resultant = {total}')
    
##############################################################################
# constants
gidMass = 0.350 #grams
fullBin = 252 # number of GID is a full bin
#theta = 2.48*np.pi/180 # position of primary link whgen full
theta = 32.345*np.pi/180 # position of primary link whgen full
FOS = 2 #factor of safety
g = 9.81
rho = 8000 # density if stainless steel AISI 316 Stainless Steel (kg/m3)
P = gidMass*fullBin*g*FOS # permitted 
pumpCapacity = 12*1000

#l1=0.58362
#L1_intial = 0.58362 # connecting link length
iter_l1 = 100
L1_list = np.linspace(0.1,0.5,num=iter_l1)


#l2 = 0.12554 # fixing pointfrom D
iter_l2 = 100
L2_list = np.linspace(0.05,0.5,num=iter_l2)

#l3 = 0.095352 # veritical displacement
iter_l3 = 100
L3_list = np.linspace(0.1,0.17,num=iter_l3)

count = 0
state = True
validity = True
requiredForce = 0
alpha = 0
l1_array = list()
l2_array = list()
l3_array = list()
count_array = list()
forces = list()
alpha_ang = list()
data = list()
target = list()
targetForce = 0
print('Running analysis..!\n')
for l1 in tqdm (L1_list, desc="Loading..."):
    for l3 in L3_list:
        for l2 in L2_list:
            ##############################################################################
            ''' 
            Geometric propteries:
            link length based on material volume and density
            '''
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
            # secondary variables
            l4 = (l_CD/2)-l2
            argument = (l3 + l2*np.sin(theta))/l1
            
            if argument <= 1:
                alpha = np.arcsin(argument)
            else:
                math.isnan(alpha) == True

            if math.isnan(alpha) == True:
                validity = False
            else:
                #print(alpha)
                ##############################################################################
                # calculation
                #Link B
                Rbx = 0
                Rby = (0.5*l_AB*(P+w_AB)-(l_CD*np.cos(theta)*(P+w_AB)))/(l_CD*np.cos(theta))
                Rcx = 0
                Rcy = w_AB + P - Rby

                # Link EB
                Rfx = 0
                Rfy = -2*Rby - w_EB
                Rex = 0
                Rey = Rfy + w_EB + Rby 

                #Link CD
                F_hp = (0.5*l_CD*(2*Rcy - Rfy + w_CD))/(-0.5*l_CD*np.sin(alpha) + (l4-0.5*l_CD)*np.tan(theta)*np.cos(alpha) + l4*np.sin(alpha))
                Rdx = - F_hp *np.cos(alpha)
                Rdy = Rcy - Rfy + w_CD - F_hp*np.sin(alpha)

                #results (outputs)
                requiredForce = round(F_hp*np.cos(alpha),2)
                angularPosition = round(theta*180/np.pi,2)
                #massOfGIDs = round((P*FOS)/g,2)
                #numberGID = round((massOfGIDs)/gidMass,0)
                targetForce = round(F_hp*np.sin(alpha),2)

            if abs(requiredForce) > abs(pumpCapacity/FOS) :
                #print(requiredForce)
                #print(count)
                pass
            else:
                newList = []
                #print(f'{count} - L1: {round(l1,4)}, L2: {round(l2,4)}, L3: {round(l3,4)}, {requiredForce}')

                l1_array.append(round(l1,4))
                l2_array.append(round(l2,4))
                l3_array.append(round(l3,4))
                forces.append(abs(requiredForce))
                count_array.append(count)
                alpha_ang.append(argument*180/np.pi)
                target.append(targetForce)
                #newList = [round(l1,4),round(l2,4),round(l3,4)]
                #print(newList)
                #var = round(abs(requiredForce))
                #print(var)
                #data[var] = newList              
                #print(data)
                count+=1

                #state = False
        # if state == False:
        #     break          

    # if state == False:
    #     break      
        #print(f'{count} - Required pump force: {requiredForce} at angle {angularPosition} with {numberGID} ({massOfGIDs} kg) GIDs')

#print(requiredForce)
print('##############################################\n')
print('Final Result\n')
print(f'{max(count_array)} different possibilities found!\n')
print(f'Range L1: ({min(l1_array)} to {max(l1_array)})')
print(f'Range L2: ({min(l2_array)} to {max(l2_array)})')
print(f'Range L3: ({min(l3_array)} to {max(l3_array)})')
print('\n')
print('##############################################\n')

var = forces.index(min(forces))
print('Best values!')
#print(forces[var])
print(f'L1 - {l1_array[var]}')
print(f'L2 - {l2_array[var]}')
print(f'L3 - {l3_array[var]}')
print(f'Alpha - {alpha_ang[var]}')
print(f'force - {target[var]}')
print('##############################################\n')
