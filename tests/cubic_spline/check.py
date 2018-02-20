import numpy as np
fname="vel_1410_dns.prof"
data=np.loadtxt(fname,skiprows=14,usecols=(0,1,2,3,4,5,6))
y=data[:,0]
yplus=data[:,1]
Uplus=data[:,2]
urmsplus=data[:,3]
vrmsplus=data[:,4]
wrmsplus=data[:,5]
uvplus=data[:,6]

from scipy.interpolate import CubicSpline
cs=CubicSpline(yplus,uvplus)
print(cs(79.0))
exit()

from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
plt.style.use('sjc')

fig=plt.figure()
ax=fig.gca()
ax.semilogx(yplus,Uplus,'o-')
figname="Uplus_yplus.png"
plt.savefig(figname)

fig=plt.figure()
ax=fig.gca()
ax.plot(y,urmsplus,'o-',label=r"$\langle u^{\prime} u^{\prime} \rangle$")
ax.plot(y,vrmsplus,'o-',label=r"$\langle v^{\prime} v^{\prime} \rangle$")
ax.plot(y,wrmsplus,'o-',label=r"$\langle w^{\prime} w^{\prime} \rangle$")
ax.plot(y,uvplus,'o-',label=r"$\langle u^{\prime} v^{\prime} \rangle$")
ax.legend()
figname="Reynolds_Stresses_ydelta.png"
plt.savefig(figname)

fig=plt.figure()
ax=fig.gca()
ax.semilogx(yplus,urmsplus,'o-',label=r"$\langle u^{\prime} u^{\prime} \rangle$")
ax.semilogx(yplus,vrmsplus,'o-',label=r"$\langle v^{\prime} v^{\prime} \rangle$")
ax.semilogx(yplus,wrmsplus,'o-',label=r"$\langle w^{\prime} w^{\prime} \rangle$")
ax.semilogx(yplus,uvplus,'o-',label=r"$\langle u^{\prime} v^{\prime} \rangle$")
ax.legend()
ax.set_xlim([1,2E3])
ax.set_ylim([-1,3])
figname="Reynolds_Stresses_yplus.png"
plt.savefig(figname)

