import numpy as np
import matplotlib.pyplot as plt


print("Plotting ...")

g_profile1=np.loadtxt('g_out_profile1',dtype=np.float64,skiprows=1)
g_profile2=np.loadtxt('g_out_profile2',dtype=np.float64,skiprows=1)

ggt_profile1=np.loadtxt('ggt_out_profile1',dtype=np.float64,skiprows=1)
ggt_profile2=np.loadtxt('ggt_out_profile2',dtype=np.float64,skiprows=1)

fig=plt.figure(figsize=(16,4.5))
plt.subplots_adjust(wspace=0.25)
g_label=[r"$g_x$",r"$g_y$", r"$g_z$"]
linecolor=['r','g','b']
for i in range(3):
    ax=plt.subplot(1,3,i+1)
    x=g_profile1[:,0]
    plt.plot(x,g_profile1[:,i+3],'r.',label='profile1')
    plt.plot(x,g_profile2[:,i+3],'b.',label='profile2')
    plt.legend()
    plt.title(g_label[i],fontsize=16)
    plt.xlabel('x (km)',fontsize=12)
    plt.ylabel(g_label[i]+' (mGal)',fontsize=12)
plt.savefig('g.jpg', dpi=600,bbox_inches='tight')    

temp=np.array([0,1,2,4,5,8])
fig=plt.figure(figsize=(16,11))
plt.subplots_adjust(wspace=0.35,hspace=0.15)#wspace: width space, hspace: height space
ggt_label=[r"$T_{xx}$",r"$T_{xy}$", r"$T_{xz}$",r"$T_{yy}$",r"$T_{yz}$",r"$T_{zz}$"]
for i in range(6):
    ax=plt.subplot(3,3,temp[i]+1)
    plt.plot(x,ggt_profile1[:,temp[i]+3],'r.',label='profile1')
    plt.plot(x,ggt_profile2[:,temp[i]+3],'b.',label='profile2')
    plt.legend()
    if i==0 or i==3 or i==5:
        plt.xlabel('x (km)',fontsize=12)
    plt.ylabel(ggt_label[i]+r' $(s^{-2})$',fontsize=12)

plt.savefig('ggt.jpg', dpi=600,bbox_inches='tight')
# plt.show()
