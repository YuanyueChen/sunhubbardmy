#!/usr/bin/dev python
#coding=utf-8
"""
Author:         Xiao-Yan Xu <wanderxu@gmail.com>
Description:

"""
import sys
import math as mh
import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{amssymb}',r'\usepackage{wasysym}'] #for \text command

fig, ( (ax1, ax3, ax5), (ax2, ax4, ax6) ) = plt.subplots(2,3,figsize=(8, 5))
plt.rc('font', size=10)

# fig1
# y-dir lines
ax1.plot([-1,-1], [-1, 2], 'k-', lw=0.5)
ax1.plot([ 0, 0], [-1, 2], 'k-', lw=0.5)
ax1.plot([ 1, 1], [-1, 2], 'k-', lw=0.5)
ax1.plot([ 2, 2], [-1, 2], 'k-', lw=0.5)

# x-dir lines
ax1.plot([-1, 2], [-1,-1], 'k-', lw=0.5)
ax1.plot([-1, 2], [ 0, 0], 'k-', lw=0.5)
ax1.plot([-1, 2], [ 1, 1], 'k-', lw=0.5)
ax1.plot([-1, 2], [ 2, 2], 'k-', lw=0.5)

#ax1.set_xlabel( r'$1/L$', fontsize=10)
#ax1.set_ylabel( r'$m$', fontsize=10)
#ax1.set_title( r'$\square$-nesting-FS', fontsize=10)
#ax1.set_ylim([0.38,1.02])
#ax1.set_xlim([6.38,7.42])
#ax1.set_xticks([6.4,6.6,6.8,7.0,7.2,7.4])
##ax1.grid(True)

# fig3
# y-dir lines
ax3.plot([-1,-1], [-1, 2], 'k-', lw=0.5)
ax3.plot([ 0, 0], [-1, 2], 'k--', lw=0.5)
ax3.plot([ 1, 1], [-1, 2], 'k-', lw=0.5)
ax3.plot([ 2, 2], [-1, 2], 'k--', lw=0.5)

# x-dir lines
ax3.plot([-1, 2], [-1,-1], 'k-', lw=0.5)
ax3.plot([-1, 2], [ 0, 0], 'k-', lw=0.5)
ax3.plot([-1, 2], [ 1, 1], 'k-', lw=0.5)
ax3.plot([-1, 2], [ 2, 2], 'k-', lw=0.5)


# fig5
# honeycomb
ax5.plot([0,0], [0,1], 'k-', lw=0.5)
ax5.plot([0,mh.sqrt(3)/2.0], [1,1.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)/2.0, mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3), mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3), mh.sqrt(3)/2], [0, -0.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)/2, 0], [-0.5,0], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3),mh.sqrt(3)], [0,1], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3),mh.sqrt(3)*3/2.0], [1,1.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2.0, 2*mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, 2*mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2], [0, -0.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2,mh.sqrt(3)], [-0.5,0], 'k-', lw=0.5)


ax5.plot([2*mh.sqrt(3),2*mh.sqrt(3)], [0,1], 'k-', lw=0.5)
ax5.plot([2*mh.sqrt(3),mh.sqrt(3)*5/2.0], [1,1.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2.0, 3*mh.sqrt(3)], [1.5, 1.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, 3*mh.sqrt(3)], [1,0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2], [0, -0.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*2], [-0.5,0], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)/2, mh.sqrt(3)/2], [3/2,5/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)/2,mh.sqrt(3)], [5/2,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3), mh.sqrt(3)*3/2], [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2, mh.sqrt(3)*3/2], [5/2,3/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2, mh.sqrt(3)], [3/2, 1], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3), mh.sqrt(3)/2], [1,3/2], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*3/2, mh.sqrt(3)*3/2], [3/2,5/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2, 2*mh.sqrt(3)], [5/2,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, mh.sqrt(3)*5/2], [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*5/2], [5/2,3/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*2], [3/2, 1], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2], [1,3/2], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*5/2, mh.sqrt(3)*5/2], [3/2,5/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*3], [5/2,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, mh.sqrt(3)*7/2], [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2, mh.sqrt(3)*7/2], [5/2,3/2], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2, mh.sqrt(3)*3], [3/2, 1], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2], [1,3/2], 'k-', lw=0.5)



ax5.plot([mh.sqrt(3),mh.sqrt(3)],           [3,4], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3),mh.sqrt(3)*3/2.0],     [4,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2.0,2*mh.sqrt(3)],   [4.5, 4.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, 2*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2],    [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2, mh.sqrt(3)],      [2.5,3], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*2,2*mh.sqrt(3)],       [3,4], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2,mh.sqrt(3)*5/2.0],   [4,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2.0, 3*mh.sqrt(3)],  [4.5, 4.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, 3*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3, mh.sqrt(3)*5/2],    [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,2*mh.sqrt(3)],     [2.5,3], 'k-', lw=0.5)


ax5.plot([3*mh.sqrt(3),3*mh.sqrt(3)],       [3,4], 'k-', lw=0.5)
ax5.plot([3*mh.sqrt(3),mh.sqrt(3)*7/2.0],   [4,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2.0, 4*mh.sqrt(3)],  [4.5, 4.0], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*4, 4*mh.sqrt(3)],      [4,3], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2],    [3, 2.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*3],     [2.5,3], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*3/2,mh.sqrt(3)*3/2],   [4.5,5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3/2,mh.sqrt(3)*2],     [5.5,6], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2.0,mh.sqrt(3)*5/2],   [6, 5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*5/2],   [5.5,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*2],     [4.5, 4], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*2, mh.sqrt(3)*3/2],    [4,4.5], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*5/2],   [4.5,5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*5/2,mh.sqrt(3)*3],     [5.5,6], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3,  mh.sqrt(3)*7/2],   [6, 5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*7/2],   [5.5,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*3],     [4.5,4], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*3,  mh.sqrt(3)*5/2],   [4,4.5], 'k-', lw=0.5)


ax5.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*7/2],   [4.5,5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*7/2,mh.sqrt(3)*4],     [5.5,6], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*4,  mh.sqrt(3)*9/2],   [6,5.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*9/2,mh.sqrt(3)*9/2],   [5.5,4.5], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*9/2,mh.sqrt(3)*4],     [4.5,4], 'k-', lw=0.5)
ax5.plot([mh.sqrt(3)*4, mh.sqrt(3)*7/2],    [4,4.5], 'k-', lw=0.5)

plt.tight_layout()
plt.savefig("fig_lattice.pdf", dpi=300)
