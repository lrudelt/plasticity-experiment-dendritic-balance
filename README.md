# plasticity-experiment-dendritic-balance
Code to reproduce the experiment in Hayama et al. 2013 [1] using the model of voltage-dependent plasticity with dendritic balance from Mikulasch et al. 2021 [3]. Figure from [2].
![](https://pad.gwdg.de/uploads/c4ea241992182d482dec04a57.png)

The model uses the voltage-dependent plasticity rule for input synapses from [3], which implements gradient ascent on the prediction error. The prediction error is encoded by the membrane voltage $u(t)$ close to the dendritic spine of the synapse. During plasticity, the synaptic weight $F$ slowly follows the decoding weight $D$, which on faster timescales is adapted according to
$$\Delta D = \eta_D z^{\text{BAP}}(t)\frac{u(t)}{F}.$$ Here, the local membrane voltage $u(t)$ is a linear sum
$$u(t) = F x^{\text{EPSP}}(t) + W z^{\text{IPSP}}(t) - D F z^{\text{BAP}}(t)$$ of the synaptic input EPSP ($F x^{\text{EPSP}}$), inhibitory input IPSP ($W z^{\text{IPSP}}$) and a backpropagating action potential (BAP) ($DFz^{\text{BAP}}$).
The time traces are modeled as an exponential decay $$\begin{aligned}x^{\text{EPSP}}(t) &= x^{\text{EPSP}}_{\text{max}} \exp(-(t-t^{sp})/\tau^{\text{EPSP}}),\\
z^{\text{IPSP}}(t) &= z^{\text{IPSP}}_{\text{max}} \exp(-(t-t^{sp})/\tau^{\text{IPSP}}),\\ z^{\text{BAP}}(t) &= z^{\text{BAP}}_{\text{max}} \exp(-(t-t^{sp})/\tau^{\text{BAP}}).\end{aligned}$$ To match the biophysical properties of EPSP, IPSP and BAP, we note that both the maximum value of the time trace (e.g. $x^{\text{EPSP}}_{\text{max}}$), as well as the synaptic weight in the model (e.g. $F$) can be adapted. Here, we chose to set the synaptic weight values to one, i.e. $D=F=1$, and $W = -1 $. The maximum EPSP and IPSP values were then chosen to agree with the experimentally measured values (Supplementary Fig 1 in [1]):
$$\begin{aligned}
x^{\text{EPSP}}_{\text{max}} \approx 70 \,\text{pA},\\
z^{\text{IPSP}}_{\text{max}} \approx 45 \,\text{pA}.
\end{aligned}$$ Moreover, we chose the following values for the timescales of the PSPs
$$\begin{aligned}
\tau^{\text{EPSP}} = 10 \,\text{ms},\\
\tau^{\text{BAP}} = 5 \,\text{ms}.
\end{aligned}$$

Values for the spine volume changes are:
LTP: $\Delta V_H \approx 63 \pm 21 \%$
LTD: $\Delta V_H \approx 2.9 \pm 3.7 \%$
LTP: $\Delta V_H \approx - 38.0 \pm 5.2 \%$

TODO:
- [x] Mention the time trace variables
- [x] Mention why we set D and F to 1, and rather absorb experimental values in the traces
- [ ] Give parameters that were estimated from experiment, and experimental values for the values for the spine volume change
- [ ] Mention code requirements (see other repos, e.g. bei Jonas and Sebastian)

<!-- Since we do not know $F$, but only $\text{EPSP}(t) = Fx(t)$, and also do not know the strength of the BAP $z_{\text{max}}$, we set $F=D=1$ and adjusted the learning rate $\eta_D$ and the strength of the BAP $z_{\text{max}}$. The values for the EPSP and IPSPs were set to the experimental values (see below). -->

<!-- On slow timescales, then the EPSP and BAP are scaled up: -->
<!-- $$\text{BAP}  \rightarrow D^2 \cdot \text{BAP}$$ and $$\text{EPSP} \rightarrow D \cdot \text{EPSP} ,$$ leading again to a value of $D = 1$. Therefore, the synaptic weight change $\Delta F = \Delta D = D$. To achieve this with continually updated weights one could e.g. introduce integration variables for the BAP and EPSP that set the goal for the plasticity, and are triggered only when plasticity is induced, e.g. when Calcium concentration is high: -->
<!-- $$\begin{aligned} \Delta I_{\text{BAP}} &= \eta_I [\text{Ca}](t) ( D^2 \cdot \text{BAP} - I_{\text{BAP}}) \\ \Delta I_{\text{EPSP}} &= \eta_I [\text{Ca}](t) ( D \cdot \text{EPSP} - I_{\text{EPSP}})\end{aligned}$$ and $$\begin{aligned} \Delta \text{BAP} &\propto  I_{\text{BAP}} - \text{BAP} \\ \Delta \text{EPSP} &\propto  I_{\text{EPSP}} - \text{EPSP}\end{aligned}$$ -->

[1] T. Hayama et al., “GABA promotes the competitive selection of dendritic spines by controlling local Ca2+ signaling,” Nat Neurosci, vol. 16, no. 10, pp. 1409–1416, Oct. 2013, doi: 10.1038/nn.3496.

[2] F. A. Mikulasch, L. Rudelt, M. Wibral and V. Priesemann, *in preparation*.

[3] F. A. Mikulasch, L. Rudelt, and V. Priesemann, “Local dendritic balance enables learning of efficient representations in networks of spiking neurons,” PNAS, vol. 118, no. 50, Dec. 2021, doi: 10.1073/pnas.2021925118.
