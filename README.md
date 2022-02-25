# plasticity-experiment-dendritic-balance
This repository contains code for Python 3 to reproduce the experiment in Hayama et al. 2013 [1] using the model of voltage-dependent plasticity with dendritic balance from Mikulasch et al. 2021 [3]. Figure from [2].
![](https://pad.gwdg.de/uploads/c4ea241992182d482dec04a57.png)

### Model
The model uses the voltage-dependent plasticity rule for input synapses from [3], which implements gradient ascent on the prediction error. The prediction error is encoded by the membrane voltage $u(t)$ close to the dendritic spine of the synapse. During plasticity, the synaptic weight $F$ slowly follows the decoding weight <img src="https://render.githubusercontent.com/render/math?math=D">, which on faster timescales is adapted according to
<img src="https://render.githubusercontent.com/render/math?math=\frac{dD}{dt} = \eta_D z^{\text{BAP}}(t)\frac{u(t)}{F}."> Here, the local membrane voltage $u(t)$ is a linear sum
$$u(t) = F x^{\text{EPSP}}(t) + W z^{\text{IPSP}}(t) - D F z^{\text{BAP}}(t)$$ of the synaptic input EPSP ($F x^{\text{EPSP}}$), inhibitory input IPSP ($W z^{\text{IPSP}}$) and a backpropagating action potential (BAP) ($DFz^{\text{BAP}}$). To compute the change in decoder weight induced by the plasticity protocol, we integrated the learning rule for $\frac{dD}{dt}$ over 80 seconds with a stimulation rate of 1 Hz (see below for details).
Since the synaptic weight follows the decoding weight on slower timescales, after hours of plastic changes, the synaptic weight is $F=D$. Due to the choice of initial parameters ($F_0 = D_0 = 1$), the observed relative change in synaptic weight after the experiment will be $\Delta F/F_0 = \Delta D/D_0 = D-1$, which is plotted in panel C.

### Details
The time traces are modeled as an exponential decay $$\begin{aligned}x^{\text{EPSP}}(t) &= x^{\text{EPSP}}_{\text{max}} \exp(-(t-t^{\text{EPSP}})/\tau^{\text{EPSP}}),\\
z^{\text{IPSP}}(t) &= z^{\text{IPSP}}_{\text{max}} \exp(-(t-t^{\text{IPSP}})/\tau^{\text{IPSP}}),\\ z^{\text{BAP}}(t) &= z^{\text{BAP}}_{\text{max}} \exp(-(t-t^{\text{BAP}})/\tau^{\text{BAP}}).\end{aligned}$$ To match the biophysical properties of EPSP, IPSP and BAP, we note that both the maximum value of the time trace (e.g. $x^{\text{EPSP}}_{\text{max}}$), as well as the synaptic weight in the model (e.g. $F$) can be adapted. Here, we chose to set the synaptic weight values to one, i.e. $D=F=1$, and $W = -1 $. The maximum EPSP and IPSP values were then chosen to agree with the experimentally measured values relatively to each other (Supplementary Fig 1 in [1]):
$$\begin{aligned}
x^{\text{EPSP}}_{\text{max}} \approx 70 ,\\
z^{\text{IPSP}}_{\text{max}} \approx 45.
\end{aligned}$$ Moreover, we chose the following values for the timescales of the PSPs
$$\begin{aligned}
\tau^{\text{EPSP}} = 10 \,\text{ms},\\
\tau^{\text{BAP}} = 5 \,\text{ms}.
\end{aligned}$$ The relative timing of the BAP and EPSPs in the simulation was chosen similar to the experiment with $t^{\text{BAP}} = 10\,\text{ms}$, and $t^{\text{EPSP}} = 10\,\text{ms}$ for the LTP protocol and  $t^{\text{EPSP}} = 25\,\text{ms}$ for the LTD protocol. For the GABA induced IPSPs, we assumed that the temporal overlap $\int_0^T dt\, z^{\text{BAP}}(t)z^{\text{IPSP}}(t)$ with the BAP is the same as for the EPSP in the LTP protocol.
Moreover, to match the experimental results for LTP and LTD protocol, we set the magnitude of the BAP to $z_{\text{max}}^{\text{BAP}} =1$, and the learning rate $\eta_D = 3.4\cdot 10^{-5}$. Note that all units for PSPs and voltages are dimensionless, because the resistance of the membrane is not known, and plasticity changes might further be mediated via the calcium concentration, contributing yet another unknown factor.

### Experimentally measured induced spine volume changes
The experimentally measured changes in spine volume $V_H$ from [1] are
LTP: $\Delta V_H \approx 63 \pm 21 \%$
LTD: $\Delta V_H \approx 2.9 \pm 3.7 \%$
LTP: $\Delta V_H \approx - 38.0 \pm 5.2 \%$

[1] T. Hayama et al., “GABA promotes the competitive selection of dendritic spines by controlling local Ca2+ signaling,” Nat Neurosci, vol. 16, no. 10, pp. 1409–1416, Oct. 2013, doi: 10.1038/nn.3496.

[2] F. A. Mikulasch, L. Rudelt, M. Wibral and V. Priesemann, *in preparation*.

[3] F. A. Mikulasch, L. Rudelt, and V. Priesemann, “Local dendritic balance enables learning of efficient representations in networks of spiking neurons,” PNAS, vol. 118, no. 50, Dec. 2021, doi: 10.1073/pnas.2021925118.
