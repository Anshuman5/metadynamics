# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:53:46 2022

@author: akumar
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as conv
import scipy.ndimage.filters as filters
from scipy.interpolate import interp2d
savefig = True #True to save the figures
countourlabel=False
dpi = 50
figno = '1'
ylabel_fsize = 46
xlabel_fsize = 46
ytick_labelsize = 36
xtick_labelsize = 36
colorbar_fsize = 42
pad_colorbar = 0.025
marker_size = 540
marker_color_a = 'r'
marker_color_p = '#00FFFF'
marker_color_d = 'g'
xlabel='$\phi$ (rad)'
ylabel='$\psi$ (rad)'
figname = 'adp_fes'
data_amber = np.loadtxt('G_integrated_amber.txt')
data_pbe = np.loadtxt('G_integrated_pbe0.txt')
data_dftb = np.loadtxt('fes_dftbd3.dat')
def fes(data):
    unique1, counts = np.unique(data[:,0],return_counts=True)
    n1=len(unique1)
    unique2, counts = np.unique(data[:,1],return_counts=True)
    n2=len(unique2)
    fes=np.reshape(data[:,2],(n1,n2))
    mn=np.amin(fes)
    fes=(fes-mn)    
    return fes, unique1, unique2

fes_amber, unique1a, unique2a = fes(data_amber) 
Xa,Ya = np.meshgrid(unique1a,unique2a)
Za=fes_amber.reshape((len(unique1a),len(unique1a)))

fes_pbe, unique1p, unique2p = fes(data_pbe)
Xp,Yp = np.meshgrid(unique1p,unique2p)
Zp=fes_pbe.reshape((len(unique1p),len(unique1p)))*2625.5   

fes_dftb, unique1d, unique2d = fes(data_dftb)
Xd,Yd = np.meshgrid(unique1d,unique2d)
Zd=fes_dftb.reshape((len(unique1d),len(unique1d)))

plt.rc('xtick', labelsize=xtick_labelsize)
plt.rc('ytick', labelsize=ytick_labelsize)
plt.figure(figsize=(14,9.8))
plt.gca()
plt.axis([np.amin(Xa),np.amax(Xa),np.amin(Ya),np.amax(Ya)])
plt.pcolormesh(Xa,Ya,Za.T,vmin=0,vmax=120,cmap='jet')
plt.colorbar(pad=pad_colorbar).set_label(label='Energy (kJ/mol)',size=colorbar_fsize)#,weight='bold')
plt.xticks([-3,-1.5,0,1.5,3])
plt.yticks([-3,-1.5,0,1.5,3])
c_black=np.linspace(0,100,11)
g_black=c_black+5
contoursa1=plt.contour(Xa,Ya,Za.T,c_black,colors='black',linewidths=0.75)
contoursa2=plt.contour(Xa,Ya,Za.T,g_black,colors='grey',linewidths=0.75)
if countourlabel==True:
    plt.clabel(contoursa1, inline=True, fontsize=10)
    plt.clabel(contoursa2, inline=True, fontsize=10, colors='black')

plt.scatter(-2.5,2.7,c=marker_color_a,marker='o',s=marker_size)
plt.scatter(-1.45,1.,c=marker_color_a,marker='s',s=marker_size)
#plt.scatter(-1.38,1.02,c='b',marker='s',s=marker_size)
plt.scatter(1.05,-0.85,c=marker_color_a,marker='v',s=marker_size)
plt.xlabel(xlabel, fontsize=xlabel_fsize)
plt.ylabel(ylabel, fontsize=ylabel_fsize)
plt.xlim(-3.05,3.05)
plt.ylim(-3.05,3.05)
plt.tight_layout()
if (savefig==1):
    plt.savefig(figname+'_amber_'+figno+'.png',dpi=dpi, bbox_inches='tight')


marker_color_p = '#00FFFF'
plt.figure(figsize=(14,9.8))
plt.yticks([-3,-1.5,0,1.5,3])
plt.xticks([-3,-1.5,0,1.5,3])
filtermatrix = np.ones((25,25))/np.power(25,2)
Zconv=conv.convolve2d(Zp.T,filtermatrix,boundary='wrap',mode='same')
Zconv = filters.gaussian_filter(Zp.T,10,mode='wrap')
plt.axis([np.amin(Xp),np.amax(Xp),np.amin(Yp),np.amax(Yp)])
plt.pcolormesh(Xp,Yp,Zconv,vmin=0,vmax=120,cmap='jet')
plt.colorbar(pad=pad_colorbar).set_label(label='Energy (kJ/mol)',size=colorbar_fsize)#,weight='bold')
contoursp1 = plt.contour(Xp,Yp,Zconv,c_black,colors='black',linewidths=0.75)
contoursp2 = plt.contour(Xp,Yp,Zconv,g_black,colors='grey',linewidths=0.75)
if countourlabel==True:
    plt.clabel(contoursp1, inline=True, fontsize=10)
    plt.clabel(contoursp2, inline=True, fontsize=10, colors='black')
plt.scatter(-2.8,2.9,c=marker_color_p,marker='o',s=marker_size)
plt.scatter(-1.5,1.5,c=marker_color_p,marker='s',s=marker_size)
plt.scatter(1.25,-0.95,c=marker_color_p,marker='v',s=marker_size)
plt.xlabel(xlabel, fontsize=xlabel_fsize)
plt.ylabel(ylabel, fontsize=ylabel_fsize)
plt.xlim(-3.05,3.05)
plt.ylim(-3.05,3.05)
plt.tight_layout()
if (savefig==1):
    plt.savefig(figname+'_pbe_'+figno+'.png',dpi=dpi, bbox_inches='tight')



plt.figure(figsize=(14,9.8))
plt.yticks([-3,-1.5,0,1.5,3])
plt.xticks([-3,-1.5,0,1.5,3])
f = interp2d(Xd, Yd, Zd, kind='cubic')
xnew = np.linspace(-np.pi, np.pi, 510)
ynew = np.linspace(-np.pi, np.pi, 510)
data1 = f(xnew,ynew)
Xn, Yn = np.meshgrid(xnew, ynew)
plt.pcolormesh(Xn, Yn, data1, vmin=0,vmax=120,cmap='jet')#cmap='RdBu')
filtermatrixd = np.ones((25,25))/np.power(25,2)
Zconvd=conv.convolve2d(Zd.T,filtermatrixd,boundary='wrap',mode='same')
Zconvd= filters.gaussian_filter(Zd,1.2,mode='wrap')
plt.axis([np.amin(Xd),np.amax(Xd),np.amin(Yd),np.amax(Yd)])
plt.colorbar(pad=pad_colorbar).set_label(label='Energy (kJ/mol)',size=colorbar_fsize)#,weight='bold')
contoursd1 = plt.contour(Xd,Yd,Zconvd,c_black,colors='black',linewidths=0.75)
contoursd2 = plt.contour(Xd,Yd,Zconvd,g_black,colors='grey',linewidths=0.75)
if countourlabel==True:
    plt.clabel(contoursd1, inline=True, fontsize=10)
    plt.clabel(contoursd2, inline=True, fontsize=10, colors='black')
plt.scatter(-2.8,2.9,c=marker_color_p,marker='o',s=marker_size)
plt.scatter(-1.5,1.5,c=marker_color_p,marker='s',s=marker_size)
plt.scatter(1.25,-0.95,c=marker_color_p,marker='v',s=marker_size)
plt.scatter(-2.8,2.9,c=marker_color_d,marker='o',s=marker_size)
plt.scatter(-1.36,1.511,c=marker_color_d,marker='s',s=marker_size)
plt.scatter(1.25,-1.15,c=marker_color_d,marker='v',s=marker_size)
plt.scatter(-2.5,2.7,c=marker_color_a,marker='o',s=marker_size)
plt.scatter(-1.45,1.,c=marker_color_a,marker='s',s=marker_size)
#plt.scatter(-1.38,1.02,c='b',marker='s',s=marker_size)
plt.scatter(1.05,-0.85,c=marker_color_a,marker='v',s=marker_size)
plt.xlabel(xlabel, fontsize=xlabel_fsize)
plt.ylabel(ylabel, fontsize=ylabel_fsize)
plt.tight_layout()
plt.xlim(-3.05,3.05)
plt.ylim(-3.05,3.05)
plt.show()
if (savefig==True):
    plt.savefig(figname+'_dftb_'+figno+'.png',dpi=dpi, bbox_inches='tight')


