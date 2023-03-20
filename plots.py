import numpy as np
import matplotlib.pyplot as plt

#Init variables
pulsars=[]
P_list,P_dot_list,x_s_list,y_s_list,em_type_list=[],[],[],[],[]
P_str,P_dot_str,x_s_str,y_s_str,emission_type='','','','',''
P,P_dot,x_s,y_s=0,0,0,0

#read the file in the pulsars list
file= open("P_Pdot_positions.txt","r")
content= file.readlines()

for i in range(len(content)):
    pulsars+=[content[i]]

for i in range(len(pulsars)):
    for j in range(12):
        P_str+=pulsars[i][j]
        P_dot_str+=pulsars[i][j+13]
        x_s_str+=pulsars[i][j+27]
        y_s_str+=pulsars[i][j+39]
    emission_type+=pulsars[i][52]+pulsars[i][53]
    P=float(P_str)
    P_dot=float(P_dot_str)
    x_s=float(x_s_str)
    y_s=float(y_s_str)
    pulse=[P,P_dot,x_s,y_s,emission_type]
    pulsars[i]=pulse
    P_str,P_dot_str,x_s_str,y_s_str,emission_type='','','','',''

#Creation of the different useful lists
for i in range(len(pulsars)):
    P_list+=[pulsars[i][0]]
    P_dot_list+=[pulsars[i][1]]
    x_s_list+=[pulsars[i][2]]
    y_s_list+=[pulsars[i][3]]
    em_type_list+=[pulsars[i][4]]

#P-Pdot plot
plt.figure(1)
plt.scatter(P_list,P_dot_list,c='red',marker='o',s=20)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
plt.title("Spin period derivative - Spin period diagram")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot.png')

#Positions plot
#plt.figure(2)
#plt.scatter(x_s_list,y_s_list)
#plt.scatter([0],[0],c='yellow',marker='o',s=40) #position of the sun
#plt.xlim(-50,30)
#plt.ylim(-20,25)
#plt.title('Positions of the detected pulsars compared to the sun')
#plt.xlabel('x (kpc)')
#plt.ylabel('y (kpc)')
#plt.savefig('Positions_detected_pulsars.png')

