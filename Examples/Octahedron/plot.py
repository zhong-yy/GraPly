import numpy as np
import matplotlib.pyplot as plt


print("Plotting ...")

g_profile1=np.loadtxt('g_out_profile1',dtype=np.float64,skiprows=1)
g_profile2=np.loadtxt('g_out_profile2',dtype=np.float64,skiprows=1)

ggt_profile1=np.loadtxt('ggt_out_profile1',dtype=np.float64,skiprows=1)
ggt_profile2=np.loadtxt('ggt_out_profile2',dtype=np.float64,skiprows=1)

fig=plt.figure(figsize=(16,4.5))
plt.subplots_adjust(wspace=0.3)
g_label=[r"$g_x$",r"$g_y$", r"$g_z$"]
linecolor=['r','g','b']
for i in range(3):
    ax=plt.subplot(1,3,i+1)
    x=g_profile1[:,0]
    plt.plot(x,g_profile1[:,i+3],'rd',markersize=6,label='Profile 1: outside the octahedron')
    plt.plot(x,g_profile2[:,i+3],'b.',markersize=12,label='Profile 2: passing through the octahedron')
#    plt.title(g_label[i],fontsize=18)
    plt.xlabel('x (km)',fontsize=16)
    plt.ylabel(g_label[i]+' (mGal)',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

ax=plt.subplot(1,3,1)
plt.legend(bbox_to_anchor=(-0.02, 1.),ncol=2,loc='lower left',fontsize=16)


plt.savefig('g.jpg', dpi=600,bbox_inches='tight')    

temp=np.array([0,1,2,4,5,8])
fig=plt.figure(figsize=(16,11))
plt.subplots_adjust(wspace=0.35,hspace=0.25)#wspace: width space, hspace: height space
ggt_label=[r"$T_{xx}$",r"$T_{xy}$", r"$T_{xz}$",r"$T_{yy}$",r"$T_{yz}$",r"$T_{zz}$"]
for i in range(6):
    ax=plt.subplot(3,3,temp[i]+1)
    plt.plot(x,ggt_profile1[:,temp[i]+3],'rd',markersize=6,label='Profile 1: outside the octahedron')
    plt.plot(x,ggt_profile2[:,temp[i]+3],'b.',markersize=12,label='Profile 2: passing through the octahedron')

    if i==0 or i==3 or i==5:
        plt.xlabel('x (km)',fontsize=16)
    plt.ylabel(ggt_label[i]+r' $(s^{-2})$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    ax.yaxis.offsetText.set_fontsize(16)

ax=plt.subplot(3,3,5)
plt.legend(bbox_to_anchor=(-1.2, -0.5),loc='upper left',fontsize=20)
plt.savefig('ggt.jpg', dpi=600,bbox_inches='tight')
# plt.show()
