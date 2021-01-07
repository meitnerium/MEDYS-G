import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,
suffixe="_2CSF_800nm_I9_10_14_R4_0"
xlabeltime="Time (u.a.)"
ylabelchamp="Intensity (u.a.)"
ylabelpop="Population"
panelxpos=0.85
xlabelspectre="k"
ylabelspectre="Intensity"
tmin=0
tmax=3E3
popmax1=1.2E-3
popmax2=0.4E-3
FILECHAMP='champ'
REP80=''
champ80=np.loadtxt( 'champ.dat')

fig = plt.figure() 
fig.canvas.set_window_title('CSFQ')
plt.subplot(3, 1, 1)
plt.plot(champ80[:,0],champ80[:,2],"r-")
plt.xlabel(xlabeltime)
plt.ylabel(ylabelchamp)
plt.figtext(panelxpos,0.80,"(A)")



panelxpos=0.75
tmin=0
tmax=3E3
popmax1=1.E0
popmax2=1.E0
pop1=np.genfromtxt('CSF_Q.dat')
print(pop1.shape)
print(pop1.shape)
plt.subplot(3, 1, 2)
plt.plot(champ80[:,0],pop1[:,0],"b-")
plt.axis([champ80[0,0], champ80[champ80[:,0].size-1,0], np.min(pop1[:,0]), 1.])
plt.xlabel(xlabeltime)
plt.ylabel(ylabelpop)

plt.subplot(3, 1, 3)
tn=np.arange(0., pop1.size/3, 1.)
print(tn.shape)
print(pop1.shape)
plt.plot(champ80[:,0],pop1[:,1],"g--",pop1[:,2],"g-")
plt.axis([tmin, champ80[champ80[:,0].size-1,0], -0.01, np.max(pop1[:,1])])
plt.xlabel(xlabeltime)
plt.ylabel(ylabelpop)

plt.tight_layout()
fig.savefig( 'CSF_Q'+suffixe+'.png')
print('CSF_Q.png created')

CSFPop=np.genfromtxt('CSF_P.dat')
fig2 = plt.figure() 
fig2.canvas.set_window_title('CSFP')
#plt.subplot(3, 1, 1)
plt.plot(champ80[:,0],CSFPop[:,0],'b-',champ80[:,0],CSFPop[:,1],'r-',champ80[:,0],CSFPop[:,1]+CSFPop[:,0],'g-')
#plt.plot(champ80[:,0],CSFPop[:,0],"b-", CSFPop[:,1], "r-")
#plt.plot(champ80[:,0],CSFPop[:,0],"b-", CSFPop[:,1], 'r-', CSFPop[:,0]+CSFPop[:,1],'g-')
plt.axis([np.min(champ80[:,0]), np.max(champ80[:,0]), np.min(CSFPop[:,0]), np.max(CSFPop[:,0]+CSFPop[:,1])])
plt.ylabel(ylabelpop)
plt.xlabel(xlabeltime)

#plt.subplot(3, 1, 2)
#plt.plot(champ80[:,0],CSFPop[:,1],'r-',CSFPop[:,1]+CSFPop[:,0],'g-')
#plt.axis([np.min(champ80[:,0]), np.max(champ80[:,0]), np.min(CSFPop[:,1]), np.max(CSFPop[:,1])])
#plt.xlabel(xlabeltime)
#plt.ylabel(ylabelpop)

#plt.subplot(3, 1, 3)
#plt.plot(champ80[:,0],CSFPop[:,1]+CSFPop[:,0],'g-')
#plt.axis([np.min(champ80[:,0]), np.max(champ80[:,0]), np.min(CSFPop[:,1]), np.max(CSFPop[:,1])])
#plt.xlabel(xlabeltime)
#plt.ylabel(ylabelpop)

#plt.subplot(2, 2, 4)
#x2=np.loadtxt('result_800nm_without_GS/R_1.40/spectre_canal_00002.dat')
#plt.plot(pop2[:,1])
#plt.axis([tmin, tmax, 0, popmax2])
#plt.xlabel(xlabeltxt)
#plt.ylabel(ylabeltxt)
#plt.figtext(panelxpos2,panelypos2,"(D)")
#TODO DEL hspace
#plt.hspace('0')
plt.tight_layout()

CSFPFILE= 'CSF_P'+suffixe+'.png'
plt.savefig(CSFPFILE)
print(CSFPFILE + ' created')



fig4a = plt.figure() 
fig4a.canvas.set_window_title('spectre canal 1')
canal1_2d=np.loadtxt('spectre2d_canal00001.dat')
canal2_2d=np.loadtxt('spectre2d_canal00002.dat')
X=np.loadtxt('kx.dat')
dx=X[0,1]-X[0,0]
print("dx=",dx)
Y=np.loadtxt('ky.dat')
dy=Y[1,0]-Y[0,0]
print("dy=",dy)
leveltot = plt.MaxNLocator(nbins=15).tick_values(canal1_2d.min(), canal1_2d.max())
cf=plt.contourf(X,Y ,canal1_2d,levels=leveltot)
plt.axis([-2, 2, -2, 2])
fig4a.colorbar(cf)
FILESPECTRE2Dcanal1='spectre2dcanal1'+suffixe+'.png'
plt.ylabel('kz')
plt.xlabel('kx')
plt.savefig(FILESPECTRE2Dcanal1)

fig4b = plt.figure() 
fig4b.canvas.set_window_title('spectre canal 1LOG')
Z=np.log10(canal1_2d)
leveltot = plt.MaxNLocator(nbins=15).tick_values(Z.min(), Z.max())
cfb=plt.contourf(X,Y ,Z,levels=leveltot)
plt.axis([-15, 15, -15, 15])
fig4b.colorbar(cfb)
FILESPECTRE2Dcanal1LOG='spectre2dcanal1LOG'+suffixe+'.png'
plt.ylabel('kz')
plt.xlabel('kx')
plt.savefig(FILESPECTRE2Dcanal1LOG)





#fig, (ax0, ax1) = plt.subplots(nrows=2)
fig5a = plt.figure() 
fig5a.canvas.set_window_title('spectre canal 2')
leveltot = plt.MaxNLocator(nbins=15).tick_values(canal2_2d.min(), canal2_2d.max())
#cmap = plt.get_cmap('bgr')
#cmap = plt.get_cmap('PiYG')
#cf=plt.contourf(X[:-1, :-1] + dy/2.,Y[:-1, :-1] + dy/2. ,
#Z,levels=levels)
cf=plt.contourf(X,Y ,canal2_2d,levels=leveltot)
plt.axis([-2, 2, -2, 2])
fig5a.colorbar(cf)
plt.ylabel('kz')
plt.xlabel('kx')
FILESPECTRE2Dcanal2='spectre2dcanal2'+suffixe+'.png'
plt.savefig(FILESPECTRE2Dcanal2)

#fig, (ax0, ax1) = plt.subplots(nrows=2)
fig5b = plt.figure() 
fig5b.canvas.set_window_title('spectre canal 2LOG')
#cmap = plt.get_cmap('bgr')
#cmap = plt.get_cmap('PiYG')
#cf=plt.contourf(X[:-1, :-1] + dy/2.,Y[:-1, :-1] + dy/2. ,
#Z,levels=levels)
Z=np.log10(canal2_2d)
leveltot = plt.MaxNLocator(nbins=15).tick_values(Z.min(), Z.max())
cf=plt.contourf(X,Y ,Z,levels=leveltot)
plt.axis([-15, 15, -15, 15])
fig5b.colorbar(cf)
plt.ylabel('kz')
plt.xlabel('kx')
FILESPECTRE2DLOGcanal2='spectre2dcanal2LOG'+suffixe+'.png'
plt.savefig(FILESPECTRE2DLOGcanal2)


canaltot_2d=np.loadtxt('spectre_tot_coer.dat')
Z=np.log10(canaltot_2d)
fig7a=plt.figure()
fig7a.canvas.set_window_title('spectre total')
leveltot = plt.MaxNLocator(nbins=15).tick_values(canaltot_2d.min(), canaltot_2d.max())
#cmap = plt.get_cmap('bgr')
#cmap = plt.get_cmap('PiYG')
#cf=plt.contourf(X[:-1, :-1] + dy/2.,Y[:-1, :-1] + dy/2. ,
#Z,levels=levels)
cf=plt.contourf(X,Y ,canaltot_2d,levels=leveltot)
plt.axis([-2, 2, -2, 2])
plt.ylabel('kz')
plt.xlabel('kx')
fig7a.colorbar(cf)
FILESPECTRE2DTOT='spectre2dtot'+suffixe+'.png'
plt.savefig(FILESPECTRE2DTOT)

canaltot_2d=np.loadtxt('spectre_tot_coer.dat')
Z=np.log10(canaltot_2d)
fig7b=plt.figure()
fig7b.canvas.set_window_title('spectre total')
leveltot = plt.MaxNLocator(nbins=15).tick_values(Z.min(), Z.max())
#cmap = plt.get_cmap('bgr')
#cmap = plt.get_cmap('PiYG')
#cf=plt.contourf(X[:-1, :-1] + dy/2.,Y[:-1, :-1] + dy/2. ,
#Z,levels=levels)
cfb=plt.contourf(X,Y ,Z,levels=leveltot)
plt.axis([-15, 15, -15, 15])
plt.ylabel('kz')
plt.xlabel('kx')
fig7b.colorbar(cfb)
FILESPECTRE2DTOT='spectre2dtotLOG'+suffixe+'.png'
plt.savefig(FILESPECTRE2DTOT)




plt.show()
plt.plot()


