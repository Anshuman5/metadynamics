# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 12:14:00 2022

@author: akumar
"""
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import scipy.signal as conv
import scipy.ndimage.filters as filters
savefig = 0 #True to save the figures
#savefig = True
countourlabel=True
dpi = 50
figno = '5'
ylabel_fsize = 46
xlabel_fsize = 46
ytick_labelsize = 36
xtick_labelsize = 36
colorbar_fsize = 40
colorbar_labelsize = 27
contour_fsize = 15
pad_colorbar = 0.025
marker_size = 200
annotation_fontsize = 36
annotation_color = 'lime'
xlabel = '$\phi$ (rad)'
ylabel = '$\psi$ (rad)'
fontweight = 'normal'
# all_labels_a = [r'A$^{a}$', 'B$^{a}$']
# all_labels_d = [r'A$^{d}$', 'B$^{d}$']
all_labels_a = [r'A', 'B']
all_labels_d = [r'C', 'D']
marker_color_a = '#FF0000'
marker_color_p = '#00FFFF'
color_y = '#FFFF00'
color_g = '#FF7F00'
vmax = 60

#markers for minima
phi11 = [-2.242, -0.42, -2.20, -0.24]
psi11 = [-1.964, 0.746, -0.94, 0.76]
annotate_xy_A= [-2.39, -1.05, -2.46, -1.99]
annotate_xy_B = [-0.43, 0.72, -0.57, 0.75]
data_file = np.loadtxt('fes_remdesivir_amber.dat')
x = data_file[:,0]
y = data_file[:,1]
z = data_file[:,2]
ycli = np.linspace(3.14,-3.14)
xcli = np.linspace(-3.14, 3.14)

X, Y = np.meshgrid(xcli, ycli)
Z = griddata((x, y), z, (X, Y), method='cubic')

plt.figure(figsize=(14,9.8))
plt.gca()
plt.rc('xtick', labelsize=xtick_labelsize)
plt.rc('ytick', labelsize=ytick_labelsize)
# contours = plt.contour(X, Y, Z, 10, alpha=1, linewidths=2, linestyles='solid', extend='neither', cmap = 'Greys')
#c_black=np.linspace(0,90,11)
#g_black=c_black+5
contours1 = plt.contour(X, Y, Z, colors='k', alpha=0.5, linewidths=0.75, linestyles='solid', extend='neither')
#contours2 = plt.contour(X, Y, Z, g_black, colors='grey', alpha=0.5, linewidths=0.75, linestyles='solid', extend='neither')
if countourlabel==True:
    plt.clabel(contours1, inline=True, fontsize=contour_fsize, fmt='%1.2f')
 #   plt.clabel(contours2, inline=True, fontsize=contour_fsize)
plt.imshow(Z, extent=[ -3.14, 3.14, -3.14, 3.14], cmap = 'jet', interpolation='bicubic', alpha=1, aspect='auto', vmin=0, vmax=vmax)
plt.xlabel(xlabel, fontsize=xlabel_fsize, fontname = 'Arial', fontweight=fontweight)
plt.ylabel(ylabel, fontsize=ylabel_fsize, fontname = 'Arial', fontweight=fontweight)
cbar1 = plt.colorbar(orientation="vertical", shrink=1, pad=0.025)
cbar1.ax.tick_params(labelsize=colorbar_labelsize)
cbar1.set_ticks(list(np.arange(0,vmax+10,10)))
cbar1.set_label('Energy (kJ/mol)', fontname = 'Arial', fontsize = colorbar_fsize, fontweight=fontweight)
plt.yticks([-3.0, -1.5, 0, 1.5, 3.0], size=ytick_labelsize)
plt.xticks([-3.0, -1.5, 0, 1.5, 3.0], size=xtick_labelsize)

plt.ylim(-3.12, 3.12)
plt.xlim(-3.12, 3.12)
#plt.title('Remdesivir (Amber)', fontsize=20)
pathway_file = np.loadtxt("transition_path_remde_amber.txt")
phi31 = pathway_file[:,0]
psi31 = pathway_file[:,1]
E11 = pathway_file[:,2]
#plt.scatter(phi31, psi31, E11, marker='o', c=color_y)
plt.plot(phi31, psi31, marker='o', c=color_y, lw=2)
plt.annotate(all_labels_a[0], xy=(annotate_xy_A[0],annotate_xy_A[1]), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color)
plt.annotate(all_labels_a[1], xy=(annotate_xy_B[0],annotate_xy_B[1]), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color)
plt.scatter(phi11[2],psi11[2], c=marker_color_a, marker='o', s=marker_size)
plt.scatter(phi11[3],psi11[3], c=marker_color_p, marker='s', s=marker_size)

plt.tight_layout()


if (savefig==1):
    plt.savefig('fes_remdesivir_amber'+figno+'.png',dpi=dpi, bbox_inches='tight')


#plot for Remdesivir using DFTB 10ns FES
data = np.loadtxt('fes_10ns_remde_dftb.dat')
unique1, counts = np.unique(data[:,0],return_counts=True)
n1=len(unique1)
unique2, counts = np.unique(data[:,1],return_counts=True)
n2=len(unique2)
fes=np.reshape(data[:,2],(n1,n2))

X,Y = np.meshgrid(unique1,unique2)
Z=fes.reshape((len(unique1),len(unique1)))
c_black=np.linspace(0,90,11)
g_black=c_black+5

plt.figure(figsize=(14,9.8))
plt.gca()
plt.yticks([-3.0, -1.5, 0, 1.5, 3.0], size=ytick_labelsize)
plt.xticks([-3.0, -1.5, 0, 1.5, 3.0], size=xtick_labelsize)
plt.rc('xtick', labelsize=xtick_labelsize)
plt.rc('ytick', labelsize=ytick_labelsize)
plt.xlabel(xlabel, fontsize=xlabel_fsize, fontname = 'Arial', fontweight=fontweight)
plt.ylabel(ylabel, fontsize=ylabel_fsize, fontname = 'Arial', fontweight=fontweight)
f = interp2d(X, Y, Z, kind='cubic')
xnew = np.linspace(-np.pi, np.pi, 510)
ynew = np.linspace(-np.pi, np.pi, 510)
data1 = f(xnew,ynew)
Xn, Yn = np.meshgrid(xnew, ynew)
plt.axis([np.amin(X),np.amax(X),np.amin(Y),np.amax(Y)])
plt.pcolormesh(Xn,Yn,data1,vmin=0,vmax=vmax,cmap='jet')
#plt.colorbar(label='kJ/mol')
cbar2 = plt.colorbar(orientation="vertical", shrink=1, pad=0.025)
cbar2.ax.tick_params(labelsize=colorbar_labelsize)
cbar2.set_ticks(list(np.arange(0,vmax+10,10)))
cbar2.set_label('Energy (kJ/mol)', fontname = 'Arial', fontsize = colorbar_fsize, fontweight=fontweight)

filtermatrix = np.ones((25,25))/np.power(25,2)
Zconv = conv.convolve2d(data1.T,filtermatrix,boundary='wrap',mode='same')
Zconv = filters.gaussian_filter(data1,1.2,mode='wrap')
contoursd1 = plt.contour(xnew,ynew,Zconv,c_black,colors='k',alpha=0.5, linewidths=0.75, linestyles='solid', extend='neither')
contoursd2 = plt.contour(xnew,ynew,Zconv,g_black,colors='k',alpha=0.5, linewidths=0.75, linestyles='solid', extend='neither')
if countourlabel==True:
    plt.clabel(contoursd1, inline=True, fontsize=contour_fsize, fmt='%1.2f')
    plt.clabel(contoursd2, inline=True, fontsize=contour_fsize, colors='black', fmt='%1.2f') 
plt.ylim(-3.12, 3.12)
plt.xlim(-3.12, 3.12)
#plt.title('Remdesivir (DFTB)', fontsize=20)
pathway_file = np.loadtxt("test_pathway_remdesivir.txt")
phi21 = pathway_file[:,0]
psi21 = pathway_file[:,1]
E11 = pathway_file[:,2]
#plt.scatter(phi21, psi21, E11, marker='o', c="magenta")
plt.plot(phi21, psi21, marker='o', c="r", lw=2)
plt.plot(phi31, psi31, marker='o', c=color_y, lw=2) #comment 
plt.tight_layout()
plt.annotate(all_labels_d[0], xy=(annotate_xy_A[2],annotate_xy_A[3]), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color)
plt.annotate(all_labels_d[1], xy=(annotate_xy_B[2],annotate_xy_B[3]), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color)
plt.annotate(all_labels_a[0], xy=(annotate_xy_A[0],annotate_xy_A[1]), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color) #
plt.annotate(all_labels_a[1], xy=(annotate_xy_B[0]+0.16,annotate_xy_B[1]+0.06), ha='center', va='bottom', fontsize=annotation_fontsize, c=annotation_color) #


plt.scatter(phi11[0],psi11[0], c=color_y, marker='o', s=marker_size) #comment
plt.scatter(phi11[1],psi11[1], c=color_g, marker='s', s=marker_size) #comment
plt.scatter(phi11[2],psi11[2], c=marker_color_a, marker='o', s=marker_size)
plt.scatter(phi11[3],psi11[3], c=marker_color_p, marker='s', s=marker_size)

if (savefig==0):
    plt.savefig('fes_remdesivir_dftb_comapre'+figno+'.png',dpi=dpi, bbox_inches='tight')