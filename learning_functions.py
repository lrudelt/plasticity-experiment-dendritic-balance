import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl

#mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
#mpl.rcParams['axes.spines.bottom'] = False

folder = ""

# dFji for very low rate

Fij0 = 0.36787091
Wjk = 0.2
tau = 10.0 # [ms]
T = 10 * tau
dt = 0.01 * tau
scaling = 0.5

txi = 20.0
tzj = 10.0
tzk = 10.0

eta = 0.01


# exponential spike kernel
spike_kernel = lambda t, t0: (t>t0) * np.exp(-(t - t0)/tau)

# computes spike overlap
def so(t1, t2, dt, T):
    res = 0
    for t in np.arange(0.0,T,dt):
        res += spike_kernel(t, t1) * spike_kernel(t, t2)
    return res


# only xi zj interaction
dw1 = lambda F, txi, tzj: so(txi, tzj, dt, T) - F * so(tzj, tzj, dt, T)

# including zk
dw2 = lambda F, W, txi, tzj, tzk: dw1(F, txi, tzj) - W * so(tzj, tzk, dt, T)

# compute results
ts = np.arange(-10,13*5-10,5)
n = len(ts)
Fs_no_inhibition = np.zeros(n)
Fs_with_inhibition = np.zeros(n)
Fs_no_inhibition[0] = Fij0
Fs_with_inhibition[0] = Fij0
for i in range(n):
    if ts[i] > 0:
        dF = dw1(Fs_no_inhibition[i-1], txi, tzj)
        Fs_no_inhibition[i] = Fs_no_inhibition[i-1] + eta * dF

        dF = dw2(Fs_with_inhibition[i-1], Wjk, txi, tzj, tzk)
        Fs_with_inhibition[i] = Fs_with_inhibition[i-1] + eta * dF
    elif i > 0:
        Fs_no_inhibition[i] = Fs_no_inhibition[i-1]
        Fs_with_inhibition[i] = Fs_with_inhibition[i-1]

print(Fs_no_inhibition)
print(Fs_with_inhibition)


# PLOT TEMPORAL EVOLUTION OF WEIGHTS
fig, axs = plt.subplots(2,1,figsize=(scaling*4,2 * scaling*3.5))

axs[0].set_title(r"no inhibition")
axs[0].set_xlabel(r"$t$ [min]")
axs[0].set_ylabel(r"$F$ %")
axs[0].set_ylim((0.0,1.2))
#axs[0].set_xticks([-3*tau,0,3*tau])
#axs[0].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[0].axvline(0,color='gray',linestyle="--")
axs[0].plot(ts, Fs_no_inhibition / Fij0, color="black")

axs[1].set_title(r"with inhibition")
axs[1].set_xlabel(r"$t$ [min]")
axs[1].set_ylabel(r"$F$ %")
axs[1].set_ylim((0.0,1.2))
#axs[0].set_xticks([-3*tau,0,3*tau])
#axs[1].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[1].axvline(0,color='gray',linestyle="--")
axs[1].plot(ts, Fs_with_inhibition / Fij0, color="black")

plt.tight_layout()
plt.savefig(folder + "lr_xzz.svg")

"""
x=np.linspace(-5*tau,5*tau,500)
y=np.linspace(-5*tau,5*tau,500)

fig, axs = plt.subplots(1,4,figsize=(4*scaling*4,scaling*3.5))
axs[0].set_ylabel(r"$\Delta F_{ji}$")

axs[0].set_title(r"no $z_k$-spike")
axs[0].set_xlabel(r"$\Delta t_j$")
axs[0].set_ylim((-0.6*tau,0.3*tau))
axs[0].set_xticks([-3*tau,0,3*tau])
axs[0].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[0].axhline(0,color='gray',linestyle="--")
axs[0].plot(x,dw2(x), color="black")

axs[1].set_xlabel(r"$\Delta t_j$")
axs[1].set_title(r"$\Delta t_k = -2\tau$")
axs[1].set_ylim((-0.6*tau,0.3*tau))
axs[1].set_xticks([-3*tau,0,3*tau])
axs[1].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[1].axhline(0,color='gray',linestyle="--")
axs[1].plot(x,dw4(x,-2*tau), color="black")

axs[2].set_xlabel(r"$\Delta t_j$")
axs[2].set_title(r"$\Delta t_k = 0$")
axs[2].set_ylim((-0.6*tau,0.3*tau))
axs[2].set_xticks([-3*tau,0,3*tau])
axs[2].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[2].axhline(0,color='gray',linestyle="--")
axs[2].plot(x,dw4(x,0), color="black")

axs[3].set_xlabel(r"$\Delta t_j$")
axs[3].set_title(r"$\Delta t_k = 2\tau$")
axs[3].set_ylim((-0.6*tau,0.3*tau))
axs[3].set_xticks([-3*tau,0,3*tau])
axs[3].set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
axs[3].axhline(0,color='gray',linestyle="--")
axs[3].plot(x,dw4(x,2*tau), color="black")

plt.tight_layout()
plt.savefig(folder + "lr_xzz.svg")

# spike time difference

mpl.rcParams['axes.spines.left'] = False

x=np.linspace(-2*tau,8*tau,500)
tj = 2 * tau
tk = (2+3.5) * tau
plt.figure(figsize=(scaling*4,scaling*4))
plt.xlabel(r"$t$")
plt.ylim((-0.5,3+2.5))
plt.xticks([0,3*tau,6*tau], [r"0",r"$3\tau$",r"$6\tau$"])
plt.yticks([], [])
plt.axvline(0,linestyle="--")
plt.axvline(tj,linestyle="--")
plt.axvline(tk,linestyle="--")
plt.plot(x,spike(x,0)+3, color="lightcoral")
plt.plot(x,spike(x,tj)+1.5, color="black")
plt.plot(x,spike(x,tk), color="cornflowerblue")

plt.annotate(
    '', xy=(0, 3.7 + 1), xycoords='data',
    xytext=(tj, 3.7 + 1), textcoords='data',
    arrowprops={'arrowstyle': "<|-", 'color': 'black'})
plt.annotate(r"$\Delta t_j$", xy=(tj/2, 3.7+1.2),
    xytext=(tj/2, 3.7+1.2),
    horizontalalignment='center', verticalalignment='bottom')

plt.annotate(
    '', xy=(0, 3.5 + 1), xycoords='data',
    xytext=(tk, 3.5 + 1), textcoords='data',
    arrowprops={'arrowstyle': "<|-", 'color': 'black'})
plt.annotate(r"$\Delta t_k$", xy=((tk-tj)/2+tj, 3.5+1.2),
    xytext=((tk-tj)/2+tj, 3.5+1.2),
    horizontalalignment='center', verticalalignment='bottom', color="cornflowerblue")

plt.tight_layout()
plt.savefig(folder + "lr_spikes.svg")


# homeostasis

mpl.rcParams['axes.spines.left'] = True

x=np.linspace(1,9,500)
rho = 5
fig, ax = plt.subplots(1,1,figsize=(scaling*4,scaling*3))
plt.xlabel("spike-rate", labelpad=15)
plt.ylabel(r"$\Delta T$")
plt.ylim((-7,7))
plt.xlim((0,10))
plt.xticks([rho,2*rho], [r"$\rho$",r"$2\rho$"])
#plt.axvline(rho,linestyle="--")
plt.plot(x,x-rho, color="black")
plt.scatter(rho,0,color='black')

plt.annotate(
    '', xy=(0.5, 0.5 - rho + 1), xycoords='data',
    xytext=(3.5, 3.5 - rho + 1), textcoords='data',
    arrowprops={'arrowstyle': "<|-", 'color': 'black'})

plt.annotate(
    '', xy=(8.5, 8.5 - rho + 1), xycoords='data',
    xytext=(5.5, 5.5 - rho + 1), textcoords='data',
    arrowprops={'arrowstyle': "<|-", 'color': 'black'})


# set the x-spine
ax.spines['left'].set_position('zero')

# turn off the right spine/ticks
ax.spines['right'].set_color('none')
ax.yaxis.tick_left()

# set the y-spine
ax.spines['bottom'].set_position('zero')

# turn off the top spine/ticks
ax.spines['top'].set_color('none')
ax.xaxis.tick_bottom()

plt.tight_layout()
plt.savefig(folder + "lr_homeostasis.svg")


# stdp somatic


tau = 2.0
scaling = 0.5

# only xi zj interaction
dw2 = lambda t, F: tau*(1/2*(np.exp(-(2*np.maximum(0,t) -t)/tau)) - 1/2*F)

x=np.linspace(-5*tau,5*tau,500)
y=np.linspace(-5*tau,5*tau,500)

fig, ax = plt.subplots(1,1,figsize=(4.4*scaling,scaling*3.5))
ax.set_ylabel(r"$\Delta F_{ji}$")

ax.set_title(r" ")
ax.set_xlabel(r"$\Delta t_j$")
ax.set_ylim((-0.6*tau,0.6*tau))
ax.set_xticks([-3*tau,0,3*tau])
ax.set_xticklabels([r"$-3\tau$",r"0",r"$3\tau$"])
ax.axhline(0,color='gray',linestyle="--")
vij = 0.1
ax.plot(x,dw2(x, vij), color="black",alpha=0.4)
vij = 0.5
ax.plot(x,dw2(x, vij), color="black")
vij = 0.9
ax.plot(x,dw2(x, vij), color="black",alpha=0.4)

plt.annotate(r"$F_{ji}$", xy=(5.4*tau, -0.5),
    xytext=(5.7*tau, -0.5),
    horizontalalignment='left', verticalalignment='center',color='gray')
plt.annotate(r"", xy=(5.4*tau, 0.0),
    xytext=(5.4*tau, -1.0),
    horizontalalignment='left', verticalalignment='center',
    arrowprops={'arrowstyle': '<->', 'color': 'gray'})

plt.tight_layout()
plt.savefig(folder + "lr_xz.svg")
"""
