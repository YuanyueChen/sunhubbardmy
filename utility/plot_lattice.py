#!/usr/bin/dev python
#coding=utf-8
"""
Author:         Xiao-Yan Xu <wanderxu@gmail.com>
Description:
use chi square to do fitting.

"""
import sys
import math as mh
import numpy as np
import scipy.optimize as opt
#from scipy.interpolate import spline
from scipy.interpolate import CubicSpline
from scipy.interpolate import griddata
from matplotlib import gridspec
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import collections as mc


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{amssymb}',r'\usepackage{wasysym}'] #for \text command

###### define curvefunc for curve_fit
####def curvefunc (xv, *p0 ):
####        results = p0[0] + p0[1]*xv**p0[2]
####        return results
####
####def chi_square ( xdata, ydata, ydata_sigma, p0 ):
####        popt,pcov = opt.curve_fit(curvefunc, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
####        perr = np.sqrt(np.diag(pcov))
####        rchi_sq = np.sum( ( (ydata-curvefunc(xdata, *popt ) ) / ydata_sigma )**2 ) / len(ydata)
####        return popt, perr, rchi_sq
####
###### define curvefunc2 for curve_fit
####def curvefunc2 (xv, *p0 ):
####        results = p0[0] + p0[1]*xv
####        return results
####
####def chi_square2 ( xdata, ydata, ydata_sigma, p0 ):
####        popt,pcov = opt.curve_fit(curvefunc2, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
####        perr = np.sqrt(np.diag(pcov))
####        rchi_sq = np.sum( ( (ydata-curvefunc2(xdata, *popt ) ) / ydata_sigma )**2 ) / len(ydata)
####        return popt, perr, rchi_sq
####
####
#####assert len(sys.argv) == 7, "Usage: python file.py f1.dat f2.dat f3.dat f4.dat f5.dat f6.dat"
#####print "reading file "+sys.argv[1]+" ......"

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


# 倒格矢
#ax2 = fig.add_subplot(1, 2, 1, projection='3d')




plt.tight_layout()
plt.savefig("fig_lattice.pdf", dpi=300)
