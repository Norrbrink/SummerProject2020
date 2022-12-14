# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:11:29 2020

@author: Hugo Rauch
"""

import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import imageio
import os
#%%#################################################################################

c=3e8
cenx=0
ceny=0
cenz=0
rad=2

def n(x,y,z): #inputs position, gives refractive index
    if ((x-cenx)**2+(y-ceny)**2+(z-cenz)**2)**(0.5) <= rad:
        n = 1+0.01*((x-cenx)**2+(y-ceny)**2+(z-cenz)**2)
    else:
        n=1
    return n
        
def R(x,y,z):
    R=((x-0.5)**2+(y-0.5)**2+(z-0.5)**2)**(0.5)
    return R

def pol(x,y,z): #convert to polar
    if z==0:
        th=sp.pi/2
    else:
        th=sp.arctan((sp.sqrt(x**2+y**2))/(z))
    return th

def azi(x,y,z): #azimuthal coordinate
    if x==0:
        p=sp.pi/2
    else:
        p=sp.arctan(y/x)
    return p

def MultipleofTen(n):
    while n>0:
        n=n-40
    if n==0:
        return 1
    return 0
#%%#################################################################################


dt=1e-10
dx=1e-15

Q_1=sp.array([])
Q_2=sp.array([])
Q_3=sp.array([])
plane_set=sp.array([[0,0]])

Num=400


for j in range(0,10):
    for k in range(0,10):
        velocity=sp.array([0,c,0])
        displacement=sp.array([0.2*k,-4,0.2*j])      
        logic=0
        for i in range(0,Num):            
            phi=azi(*velocity)
            theta=pol(*velocity)

            velocity_phi_unit_component=sp.array([-sp.sin(phi),sp.cos(phi),0])
            velocity_theta_unit_component=sp.array([sp.cos(theta)*sp.cos(phi),sp.cos(theta)*sp.sin(phi),-sp.sin(theta)])

            theta=theta+((n(*(velocity_theta_unit_component*dx+displacement))-n(*displacement))/dx)*((3*10**8)/(n(*displacement)))*dt
            phi=phi+((n(*(velocity_phi_unit_component*dx+displacement))-n(*displacement))/dx)*((3*10**8)/(n(*displacement)))*dt

            velocity_new=sp.array([c/(n(*displacement))*sp.sin(theta)*sp.cos(phi),c/(n(*displacement))*sp.sin(theta)*sp.sin(phi),c/(n(*displacement))*sp.cos(theta)])
   
            displacement=velocity*dt+displacement
            
            velocity=velocity_new
    
            s_1=displacement[0]
            s_2=displacement[1]
            s_3=displacement[2]
        
            Q_1=np.append(Q_1,s_1)
            Q_2=np.append(Q_2,s_2)
            Q_3=np.append(Q_3,s_3)
        
            if s_2>=5:
                logic=logic+1
            
            if logic==1:
                pair=sp.array([[s_1,s_3]])
                plane_set=np.append(plane_set,pair,axis=0)
plane_set=np.delete(plane_set,0,axis=0)

#%%######################################################################################################
#DIAGRAM STUFF
#Data for three-dimensional scattered points

xdatamatrix=[]
for i in range(0,Num):
    xdatai=[]
    for j in range(0,100):
        xdatai.append(Q_1[Num*j+i])
    xdatamatrix.insert(i, xdatai)
ydatamatrix=[]
for i in range(0,Num):
    ydatai=[]
    for j in range(0,100):
        ydatai.append(Q_2[Num*j+i])
    ydatamatrix.insert(i, ydatai)
zdatamatrix=[]
for i in range(0,Num):
    zdatai=[]
    for j in range(0,100):
        zdatai.append(Q_3[Num*j+i])
    zdatamatrix.insert(i, zdatai)
#%%

xdata = Q_1
ydata = Q_2
zdata = Q_3
minmaxdata=sp.array([min(xdata),min(ydata),min(zdata),max(xdata),max(ydata),max(zdata)])

#%%
u=np.linspace(0,2*sp.pi,100)
v=np.linspace(0,2*sp.pi,100)

x=2*np.outer(np.cos(u),np.sin(v))
y=2*np.outer(np.sin(u),np.sin(v))
z=2*np.outer(np.ones(np.size(u)),np.cos(v))
#%%######################################################################################################
#Animation
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111,projection='3d')
outdir = r"C:\Users\Alexander Norrbrink\Programming\SummerProject\SPA"


for m in range(0,Num):  
    ax.clear()
    ax = plt.gca()
    ax.set_xlim([min(minmaxdata),max(minmaxdata)])
    ax.set_ylim([min(minmaxdata),max(minmaxdata)])
    ax.set_zlim([min(minmaxdata),max(minmaxdata)])
    #ax.view_init(theta, phi) #set point of view of graph
    ax.plot_surface(x, y, z, linewidth=0.0, cstride=4, rstride=4, alpha=0.3)
    ax.scatter3D(xdatamatrix[m],ydatamatrix[m],zdatamatrix[m], c=zdatamatrix[m], cmap='summer')
    fig.savefig(outdir+"Raytrace"+"{:03d}".format(m)+".png", dpi=100) #format    
                                
images = []
indirgif= r'C:\Users\Alexander Norrbrink\Programming\SummerProject'
finoutdirgif=r'C:\Users\Alexander Norrbrink\Programming\SummerProject\RayTrace.gif'
for filenames in os.listdir(indirgif):
    if filenames.endswith('.png'):
        outdirgif=os.path.join(indirgif, filenames)
        images.append(imageio.imread(outdirgif))
imageio.mimsave(finoutdirgif, images, fps=60)
#HISTOGRAM STUFF
                
#plane_x=plane_set[0:101,0]
#plane_z=plane_set[0:101,1]

#plt.scatter(plane_x,plane_z)
    
    
    