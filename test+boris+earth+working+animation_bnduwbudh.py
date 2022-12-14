# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 22:31:48 2020

@author: dae18
"""


class Electron:
    def __init__(self,init,method):
        self.method = method
        
        #Store initial conditions in form [x0,y0,z0,vx0,vy0,vz0]
        self.ke = []
        self.S0 = init
        
        #Define properties of particle
        self.m = 1.67e-27
        self.q = -1.6e-19
        self.mu = 8e22*1e-7
        
        #Define properties of field
        self.E = sp.array([0,0,0])
        
        
        #Get trajectories
        self.r,self.v = self.traj()
        
        
        
    def a(self,S):
        a = (self.q/self.m)*(self.E + sp.cross(S[3:],self.dipole(S[:3])))
        return sp.hstack((S[3:],a))
    
    
    def gyro_period(self,r):
        return (2*sp.pi*self.m)/(self.q*sp.sqrt(sp.dot(self.dipole(r),self.dipole(r))))
    
    def dipole(self,r):
        Bx=3*r[0]*self.mu*r[2]/sp.sqrt(sp.dot(r,r))**5
        By=3*r[1]*self.mu*r[2]/sp.sqrt(sp.dot(r,r))**5
        Bz=3*(r[2]**2)*self.mu/sp.sqrt(sp.dot(r,r))**5 - self.mu/sp.sqrt(sp.dot(r,r))**3
    
        return sp.array([Bx,By,Bz])
    
    
    def traj(self):
        
        vAc = 3e-4
        period = 20000
        v = self.S0[3:]
        x = self.S0[:3]
        R = sp.zeros((period,3))
        V = sp.zeros((period,3))
        dt = 1e-12
        j=0
        fig=plt.figure(figsize=(6,6))
        ax=fig.add_subplot(111,projection='3d')
        
        for instant in range(0,period):
            x += 0.5*v*dt   
            t = self.q / self.m * self.dipole(x) * 0.5 * dt
            s = 2 * t / (1+sp.dot(t,t))
            v_minus = v + self.q / (self.m * vAc) * self.E * 0.5 * dt
            v_prime = v_minus + sp.cross(v_minus,t)
            v_plus = v_minus + sp.cross(v_prime,s)
            v = v_plus + self.q / (self.m * vAc) * self.E * 0.5 * dt
            x += 0.5* v * dt
               
            R[instant] = sp.array(x)
            V[instant] = sp.array(v)
            
                
            dt = self.gyro_period(x)/20
            
        #[[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]] -> [[x1,x2,x3],[y1,y2,y3],...,[z1,z2,z3]]  
        R = sp.column_stack(R)
        Rminmax=sp.array([min(R[0]),min(R[1]),min(R[2]),max(R[0]),max(R[1]),max(R[2])])
        k=0
        r = sp.column_stack(R)
        for j in range(0,len(r)): 
            if self.MultipleofTen(j)==1: 
                ax.clear()
                ax = plt.gca()
                ax.set_xlim([min(Rminmax),max(Rminmax)])
                ax.set_ylim([min(Rminmax),max(Rminmax)])
                ax.set_zlim([min(Rminmax),max(Rminmax)])
                ax.plot_surface(earthx, earthy, earthz, linewidth=0.0, cstride=4, rstride=4) #only for dipole field
                ax.plot(R[0][:j],R[1][:j],R[2][:j], c='b')
                ax.scatter(r[j][0],r[j][1],r[j][2], c='r')
                fig.savefig(outdir+"BorisInt"+"{:03d}".format(k)+".png", dpi=90) #format
                k+=1   
            
        #[[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]] -> [[x1,x2,x3],[y1,y2,y3],...,[z1,z2,z3]]  
        return R,V
        
    
            
           
    def plot_traj(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(*self.r)#ax.plot([x_positions],[y_pos],[z_pos])
        rmax = 1.1*sp.amax(self.r)
        ax.set_xlim(-rmax,rmax)
        ax.set_ylim(-rmax,rmax)
        ax.set_zlim(-rmax,rmax)
        plt.show()
    
    def plot_ke(self):
        k = []
        for v in self.v:
            ek = 1/2 * self.m * sp.dot(v,v)
            k.append(ek)
        
        plt.plot(k[:500])
    
    def MultipleofTen(self,n):
        while n>0:
            n=n-60
        if n==0:
            return 1
        return 0
    
    
            
        
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio
import os
import warnings
outdir = r"C:\Users\Alexander Norrbrink\Programming\SummerProject\SPA"
warnings.filterwarnings("ignore", category=DeprecationWarning)

earthradius=6.4e6
earththeta=sp.linspace(0, 2 * sp.pi, 200)
earthphi=sp.linspace(0, sp.pi, 200)
earthx=sp.outer(sp.cos(earththeta), sp.sin(earthphi))*earthradius
earthy=sp.outer(sp.sin(earththeta), sp.sin(earthphi))*earthradius
earthz=sp.outer(sp.ones(sp.size(earththeta)), sp.cos(earthphi))*earthradius

initial_conditions = sp.array([2e7,0,0,0,5e7,3e7])
electron = Electron(initial_conditions,method="boris")
print(electron.v)

#electron.plot_traj()
#electron.plot_ke()

#%%

images = []
indirgif= r'C:\Users\Alexander Norrbrink\Programming\SummerProject'
finoutdirgif=r'C:\Users\Alexander Norrbrink\Programming\SummerProject\BorisInt.gif'
for filenames in os.listdir(indirgif):
    if filenames.endswith('.png'):
        outdirgif=os.path.join(indirgif, filenames)
        images.append(imageio.imread(outdirgif))
imageio.mimsave(finoutdirgif, images)

#%%
earthradius=6.4e6
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111,projection='3d')
earththeta=sp.linspace(0, 2 * sp.pi, 200)
earthphi=sp.linspace(0, sp.pi, 200)
earthx=sp.outer(sp.cos(earththeta), sp.sin(earthphi))*earthradius
earthy=sp.outer(sp.sin(earththeta), sp.sin(earthphi))*earthradius
earthz=sp.outer(sp.ones(sp.size(earththeta)), sp.cos(earthphi))*earthradius
ax.plot_surface(earthx,earthy,earthz, linewidth=0.0, cstride=4, rstride=4)
plt.show()