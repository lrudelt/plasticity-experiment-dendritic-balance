import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

import matplotlib as mpl
import pylab as plt
from matplotlib import rc

folder = ""

"""PARAMETERS AND EXPERIMENTAL MEASUREMENTS"""
# SIMULATION PARAMETERS
N_stimulation = 80
T = 1000. # in one second = 1000ms, one spike is observed, stimulation frequency is 1Hz

# Values inspired by experiment
tau_BAP = 5. # ms
tau_EPSP = 10. # ms

EPSP = 70. # pA
IPSP = -45. # pA

t_BAP = 10.0 # ms
t_EPSP_LTP = 10.0 # ms, coincides with BAP
t_EPSP_LTD = 25.0 # ms, 15ms after the BAP

dF_ratio_LTD_experiment = 2.9 #  pm 3.7 %
dF_ratio_LTP_experiment = 63.0 # pm 21 %
dF_ratio_LTD_GABA_experiment = - 38.0 # % pm 5.2 %

dF_ratio_LTD_experiment_err = 3.7
dF_ratio_LTP_experiment_err = 21.0
dF_ratio_LTD_GABA_experiment_err = 5.2

# Tuneable parameters
eta_D = 0.000034 # learning rate in pA/ms
BAP_max = 1.0 # pA

D0 = 1. # initial value of the decoding weight

# exponential spike kernel
spike_kernel = lambda t, t0, tau: (t>t0) * np.exp(-(t - t0)/tau)

def spike_overlap(t1, t2, tau1, tau2, T = 1000.): # T sets window of integration
    spike_overlap = lambda t: spike_kernel(t, t1, tau1) * spike_kernel(t, t2, tau2)
    return quad(spike_overlap, 0.0, T)[0]

# COMPUTE INTEGRATED WEIGHT CHANGE AFTER PROTOCOL

spike_overlap_BAP = spike_overlap(t_BAP, t_BAP, tau_BAP, tau_BAP, T)
spike_overlap_EPSP_LTP =  spike_overlap(t_EPSP_LTP, t_BAP, tau_EPSP, tau_BAP, T)
spike_overlap_EPSP_LTD =  spike_overlap(t_EPSP_LTD, t_BAP, tau_EPSP, tau_BAP, T)
# Overlap between IPSP and BAP same as between EPSP and BAP for LTP protocol
spike_overlap_IPSP = spike_overlap_EPSP_LTP
# spike_overlap_IPSP = spike_overlap(t_IPSP, t_BAP, tau_IPSP, tau_BAP, T)

# only EPSP interaction
dD_EPSP_LTP = lambda D: EPSP * BAP_max * spike_overlap_EPSP_LTP - D * BAP_max**2 * spike_overlap_BAP
dD_EPSP_LTD = lambda D: EPSP * BAP_max * spike_overlap_EPSP_LTD - D * BAP_max**2 * spike_overlap_BAP
# including IPSP
dD_IPSP_LTD = lambda D: dD_EPSP_LTD(D) + IPSP * BAP_max * spike_overlap_IPSP

D_EPSP_LTP = D0
D_EPSP_LTD = D0
D_IPSP_LTD = D0
for i in range(N_stimulation):
    print(D_EPSP_LTP, D_EPSP_LTD, D_IPSP_LTD)
    D_EPSP_LTP += eta_D * dD_EPSP_LTP(D_EPSP_LTP)
    D_EPSP_LTD += eta_D * dD_EPSP_LTD(D_EPSP_LTD)
    D_IPSP_LTD += eta_D * dD_IPSP_LTD(D_IPSP_LTD)

# compute relative weight change, where the final F is equal to D afer the plasticity protocol
dF_ratio_LTP = (D_EPSP_LTP-D0)/D0 * 100 # in %
dF_ratio_LTD = (D_EPSP_LTD-D0)/D0 * 100 # in %
dF_ratio_LTD_GABA = (D_IPSP_LTD-D0)/D0 * 100 # in %

print(dF_ratio_LTP, dF_ratio_LTD, dF_ratio_LTD_GABA)

"""PLOTS"""
scaling = 0.5 # Can be used to scale the figure size
fig, axs = plt.subplots(1,2,figsize=(scaling*7, scaling*2.8))
# PLOTTING PARAMETERS
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = False

rc('text', usetex=True)
mpl.rcParams["errorbar.capsize"] = 0.0

# Induced spine volume change
axs[0].set_ylabel(r"\begin{center}induced spine volume\\ change $\Delta V_H$ (\%)\end{center}")
axs[0].set_ylim((-50 ,100))
axs[0].set_xlim((0,4))
axs[0].set_xticks(())
axs[0].spines['left'].set_bounds((-50, 100))
axs[0].bar(1, dF_ratio_LTP_experiment, 0.6, yerr= dF_ratio_LTP_experiment_err, color = '0.6')
axs[0].bar(2, dF_ratio_LTD_experiment, 0.6, yerr= dF_ratio_LTD_experiment_err, color = '0.6')
axs[0].bar(3, dF_ratio_LTD_GABA_experiment, 0.6, yerr= dF_ratio_LTD_GABA_experiment_err, color = '0.6')
axs[0].plot((0,4), (0,0), lw = 0.5, color = "0.0")

# Induced synaptic weight change
axs[1].set_ylabel(r"\begin{center}induced synaptic weight\\ change $\Delta F$ (\%)\end{center}")
axs[1].set_ylim((-50 ,100))
axs[1].set_xlim((0,4))
axs[1].set_xticks(())
axs[1].spines['left'].set_bounds((-50, 100))
axs[1].bar(1, dF_ratio_LTP, 0.6, color = 'darkgreen')
axs[1].bar(2, dF_ratio_LTD, 0.6, color = 'darkgreen')
axs[1].bar(3, dF_ratio_LTD_GABA, 0.6, color = 'darkgreen')
axs[1].plot((-0.1,4), (0,0), lw = 0.5, color = "0.0")

plt.tight_layout()
plt.savefig(folder + "plasticity_experiment.svg")
