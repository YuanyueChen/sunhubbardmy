import matplotlib.pyplot as plt
import math as mh
import matplotlib
from matplotlib import cm
import numpy as np

matplotlib.rcParams.update({'font.size': 6})
matplotlib.rcParams.update({'axes.labelpad': -16})
plt.xticks([])
plt.yticks([])
#square
fig = plt.figure(figsize=plt.figaspect(2/3))


ax = fig.add_subplot(2, 3, (1,1))
ax.set_title('$(a)$', loc = 'left')
ax.set_frame_on(False)
ax.plot([-1,-1], [-1, 2], 'k-', lw=0.5)
ax.plot([ 0, 0], [-1, 2], 'k-', lw=0.5)
ax.plot([ 1, 1], [-1, 2], 'k-', lw=0.5)
ax.plot([ 2, 2], [-1, 2], 'k-', lw=0.5)

ax.plot([-1, 2], [-1,-1], 'k-', lw=0.5)
ax.plot([-1, 2], [ 0, 0], 'k-', lw=0.5)
ax.plot([-1, 2], [ 1, 1], 'k-', lw=0.5)
ax.plot([-1, 2], [ 2, 2], 'k-', lw=0.5)
ax.scatter(x = 0, y = 0, c = 'k')
ax.text(0.2, 0, 'A')
ax.set_xticks([])
ax.set_yticks([])
#ax.set_xlabel(r'$x$', fontsize = 10)
#ax.set_ylabel(r'$y$', fontsize = 10)




ax = fig.add_subplot(2, 3, (2,2))
ax.set_title('$(b)$', loc = 'left')
ax.set_frame_on(False)
ax.plot([-1,-1], [-1, 2], 'k-', lw=0.5)
ax.plot([ 0, 0], [-1, 2], 'k--', lw=0.5)
ax.plot([ 1, 1], [-1, 2], 'k-', lw=0.5)
ax.plot([ 2, 2], [-1, 2], 'k--', lw=0.5)

ax.plot([-1, 2], [-1,-1], 'k-', lw=0.5)
ax.plot([-1, 2], [ 0, 0], 'k-', lw=0.5)
ax.plot([-1, 2], [ 1, 1], 'k-', lw=0.5)
ax.plot([-1, 2], [ 2, 2], 'k-', lw=0.5)
ax.scatter(x = 0, y = 0, c = 'k')
ax.scatter(x = 1, y = 0, c = 'r')
ax.text(0.2,0,'A' )
ax.text(1.2,0,'B' )
ax.set_xticks([])
ax.set_yticks([])
#ax.set_xlabel(r'$x$', fontsize = 10)
#ax.set_ylabel(r'$y$', fontsize = 10)



ax = fig.add_subplot(2, 3, (3,3))
ax.set_frame_on(False)
ax.set_title('$(c)$', loc = 'left')
#ax.set_xlabel(r'$x$', fontsize = 10)
#ax.set_ylabel(r'$y$', fontsize = 10)
###########3####
ax.plot([0,0], [0,1], 'k-', lw=0.5)
ax.plot([0,mh.sqrt(3)/2.0], [1,1.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2.0, mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3), mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3), mh.sqrt(3)/2], [0, -0.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2, 0], [-0.5,0], 'k-', lw=0.5)


ax.plot([mh.sqrt(3),mh.sqrt(3)], [0,1], 'k-', lw=0.5)
ax.plot([mh.sqrt(3),mh.sqrt(3)*3/2.0], [1,1.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2.0, 2*mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, 2*mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2], [0, -0.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2,mh.sqrt(3)], [-0.5,0], 'k-', lw=0.5)


ax.plot([2*mh.sqrt(3),2*mh.sqrt(3)], [0,1], 'k-', lw=0.5)
ax.plot([2*mh.sqrt(3),mh.sqrt(3)*5/2.0], [1,1.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2.0, 3*mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, 3*mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2], [0, -0.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*2], [-0.5,0], 'k-', lw=0.5)


ax.plot([3*mh.sqrt(3),3*mh.sqrt(3)], [0,1], 'k-', lw=0.5)
ax.plot([3*mh.sqrt(3),mh.sqrt(3)*7/2.0], [1,1.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2.0, 4*mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*4], [1,0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2], [0, -0.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*3], [-0.5,0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*9/2,mh.sqrt(3)*4], [-0.5,0], 'k-', lw=0.5)



#########
ax.plot([mh.sqrt(3)/2, mh.sqrt(3)/2], [3/2,5/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2,mh.sqrt(3)], [5/2,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3), mh.sqrt(3)*3/2], [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2, mh.sqrt(3)*3/2], [5/2,3/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2, mh.sqrt(3)], [3/2, 1], 'k-', lw=0.5)
ax.plot([mh.sqrt(3), mh.sqrt(3)/2], [1,3/2], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*3/2, mh.sqrt(3)*3/2], [3/2,5/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2, 2*mh.sqrt(3)], [5/2,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, mh.sqrt(3)*5/2], [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*5/2], [5/2,3/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*2], [3/2, 1], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2], [1,3/2], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*5/2], [3/2,5/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*3], [5/2,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, mh.sqrt(3)*7/2], [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2, mh.sqrt(3)*7/2], [5/2,3/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2, mh.sqrt(3)*3], [3/2, 1], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2], [1,3/2], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*7/2, mh.sqrt(3)*7/2], [3/2,5/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*4], [5/2,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*9/2], [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*9/2, mh.sqrt(3)*9/2], [5/2,3/2], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*9/2, mh.sqrt(3)*4], [3/2, 1], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2], [1,3/2], 'k-', lw=0.5)


#####################
ax.plot([0,0],  [3,4], 'k-', lw=0.5)
ax.plot([0,mh.sqrt(3)/2],     [4,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2,mh.sqrt(3)],   [4.5, 4.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3),mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3),mh.sqrt(3)/2],    [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2,0],      [2.5,3], 'k-', lw=0.5)



ax.plot([mh.sqrt(3),mh.sqrt(3)],           [3,4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3),mh.sqrt(3)*3/2.0],     [4,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2.0,2*mh.sqrt(3)],   [4.5, 4.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, 2*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2],    [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2, mh.sqrt(3)],      [2.5,3], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*2,2*mh.sqrt(3)],       [3,4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2,mh.sqrt(3)*5/2.0],   [4,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2.0, 3*mh.sqrt(3)],  [4.5, 4.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, 3*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2],    [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,2*mh.sqrt(3)],     [2.5,3], 'k-', lw=0.5)


ax.plot([3*mh.sqrt(3),3*mh.sqrt(3)],       [3,4], 'k-', lw=0.5)
ax.plot([3*mh.sqrt(3),mh.sqrt(3)*7/2.0],   [4,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2.0, 4*mh.sqrt(3)],  [4.5, 4.0], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, 4*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2],    [3, 2.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*3],     [2.5,3], 'k-', lw=0.5)

#####################
ax.plot([0,mh.sqrt(3)/2],   [6, 5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2,mh.sqrt(3)/2],   [4.5,5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)/2,mh.sqrt(3)],     [5.5,6], 'k-', lw=0.5)
ax.plot([mh.sqrt(3),mh.sqrt(3)*3/2],   [6, 5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2,mh.sqrt(3)*3/2], [5.5,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2,mh.sqrt(3)],     [4.5, 4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3), mh.sqrt(3)/2],      [4,4.5], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*3/2,mh.sqrt(3)*3/2],   [4.5,5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3/2,mh.sqrt(3)*2],     [5.5,6], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2.0,mh.sqrt(3)*5/2],   [6, 5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*5/2],   [5.5,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*2],     [4.5, 4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2],    [4,4.5], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*5/2],   [4.5,5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*3],     [5.5,6], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3,  mh.sqrt(3)*7/2],   [6, 5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*7/2],   [5.5,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*3],     [4.5,4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*3,  mh.sqrt(3)*5/2],   [4,4.5], 'k-', lw=0.5)


ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*7/2],   [4.5,5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*4],     [5.5,6], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4,  mh.sqrt(3)*9/2],   [6,5.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*9/2,mh.sqrt(3)*9/2],   [5.5,4.5], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*9/2,mh.sqrt(3)*4],     [4.5,4], 'k-', lw=0.5)
ax.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2],    [4,4.5], 'k-', lw=0.5)


ax.set_xticks([])
ax.set_yticks([])
ax.scatter(x = 0, y = 0, c = 'k')
ax.scatter(x = 0, y = 1, c = 'r')
ax.text( 0.5, 0, 'A')
ax.text( 0.5, 1, 'B')

#######################
ax = fig.add_subplot(2, 3, (4,4), projection='3d')
ax.set_title('$(d)$', loc = 'left')
X = np.arange(-1, 1, 0.02)
Y = np.arange(-1, 1, 0.02)
X, Y = np.meshgrid(X,Y)
Z = -2*(np.cos(mh.pi*X)+np.cos(mh.pi*Y)) 
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel(r'$k_x$', fontsize = 10)
ax.set_ylabel(r'$k_y$', fontsize = 10)
ax.set_zlabel(r'$k_z$', fontsize = 10)


#pi-flux
ax = fig.add_subplot(2, 3, (5,5), projection='3d')

ax.set_title('$(e)$', loc = 'left')
X = np.arange(-1, 1, 0.02)
Y = np.arange(-1, 1, 0.02)
X, Y = np.meshgrid(X,Y)
Z = np.sqrt( 4*np.cos(mh.pi*Y)*np.cos(mh.pi*Y) + (1+np.cos(2*mh.pi*X))*(1+np.cos(2*mh.pi*X)) + np.sin(2*mh.pi*X)*np.sin(2*mh.pi*X) ) 
Zp = -Z

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
surfp = ax.plot_surface(X, Y, Zp, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel(r'$k_x$', fontsize = 10)
ax.set_ylabel(r'$k_y$', fontsize = 10)
ax.set_zlabel(r'$k_z$', fontsize = 10)
#fig.colorbar(surf, shrink=0.5, aspect=10)


#honeycomb
ax = fig.add_subplot(2, 3, (6,6), projection='3d')

ax.set_title('$(f)$', loc = 'left')
X = np.arange(-1, 1, 0.01) 
Y = np.arange(-1, 1, 0.01) 
X, Y = np.meshgrid(X,Y) 
Z = np.sqrt( (1+np.cos(np.sqrt(3)/2*X*mh.pi-3/2*Y*mh.pi)+np.cos(np.sqrt(3)/2*X*mh.pi+3/2*Y*mh.pi))*(1+np.cos(np.sqrt(3)/2*X*mh.pi-3/2*Y*mh.pi)+np.cos(np.sqrt(3)/2*X*mh.pi+3/2*Y*mh.pi)) + (np.sin( np.sqrt(3)/2*X*mh.pi+3/2*Y*mh.pi )-np.sin( np.sqrt(3)/2*X*mh.pi-3/2*Y*mh.pi )) * (np.sin( np.sqrt(3)/2*X*mh.pi+3/2*Y*mh.pi )-np.sin( np.sqrt(3)/2*X*mh.pi-3/2*Y*mh.pi))) 
Zp = -Z

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
surfp =  ax.plot_surface(X, Y, Zp, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel(r'$k_x$', fontsize = 10)
ax.set_ylabel(r'$k_y$', fontsize = 10)
ax.set_zlabel(r'$k_z$', fontsize = 10)
#fig.colorbar(surf, shrink=0.5, aspect=10)



#plt.tight_layout()
plt.savefig("fig_bz.pdf", dpi=300)
