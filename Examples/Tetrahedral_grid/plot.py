import numpy as np
import matplotlib.pyplot as plt

print("Plotting ...")
g=np.loadtxt('g.dat',dtype=np.float64,skiprows=1)
ggt=np.loadtxt('ggt.dat',dtype=np.float64,skiprows=1)

x=g[:,0]
y=g[:,1]

Ndx=21
Ndy=21
X=x.reshape(Ndx,Ndy)
Y=y.reshape(Ndx,Ndy)

fig=plt.figure(figsize=(16,4))
plt.subplots_adjust(wspace=0.25)
g_label=[r"$g_x$",r"$g_y$", r"$g_z$"]
for i in range(3):
    ax=plt.subplot(1,3,i+1)
    g_field=g[:,i+3]
    g_field=g_field.reshape(Ndx,Ndy)
    plt.contourf(X,Y,g_field,cmap='RdBu_r')
    plt.title(g_label[i],fontsize=14)
    plt.xlabel('x (km)',fontsize=12)
    plt.ylabel('y (km)',fontsize=12)
    plt.gca().invert_yaxis()
    cb=plt.colorbar()
    cb.ax.set_title('mGal')
plt.savefig('g.jpg', dpi=600,bbox_inches='tight')

g_field=g[:,5]
g_field=g_field.reshape(Ndx,Ndy)
gz=g_field[(Ndx+1)//2,:]
xx=X[(Ndx+1)//2,:]
fig=plt.figure()
plt.plot(xx,gz)
plt.ylabel(r'$g_z$'+' (mGal)',fontsize=12)
plt.xlabel('x (km)',fontsize=12)
plt.savefig('gz_y=0.jpg', dpi=600,bbox_inches='tight')

temp=np.array([0,1,2,4,5,8])
fig=plt.figure(figsize=(15,11))
plt.subplots_adjust(wspace=0.25,hspace=0.25)#wspace: width space, hspace: height space
ggt_label=[r"$T_{xx}$",r"$T_{xy}$", r"$T_{xz}$",r"$T_{yy}$",r"$T_{yz}$",r"$T_{zz}$"]
for i in range(6):
    ax=plt.subplot(3,3,temp[i]+1)
    g_field=ggt[:,temp[i]+3]
    g_field=g_field.reshape(Ndx,Ndy)*1e9
    plt.contourf(X,Y,g_field,cmap='RdBu_r')
    plt.title(ggt_label[i],fontsize=14)
    if i==0 or i==3 or i==5:    
        plt.xlabel('x (km)',fontsize=12)
        plt.ylabel('y (km)',fontsize=12)
    plt.gca().invert_yaxis()
    cb=plt.colorbar()
    cb.ax.set_title(r'$E$')        
plt.savefig('ggt.jpg', dpi=600,bbox_inches='tight')
