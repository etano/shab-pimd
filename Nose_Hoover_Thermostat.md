# Nose-Hoover Thermostat #

The Nose-Hoover thermostat centers on an extended system method ([Hunenberger2005](References.md)). The basic Nose-Hoover thermostat introduces one additional degree of freedom, <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{s}'><img src='http://latex.codecogs.com/png.latex?\tilde{s}%.png' title='\tilde{s}' /></a>, such that the Lagrangian becomes

<a href='http://www.codecogs.com/eqnedit.php?latex=L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})@plus;\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})'><img src='http://latex.codecogs.com/png.latex?L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})+\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})%.png' title='L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})+\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})' /></a>.

where the ~ denotes extended system variables. Now our conserved quantity is simply the original Hamiltonian plus the energy of the added Nose-Hoover degrees of freedom. It follows that the equations of motion are

<a href='http://www.codecogs.com/eqnedit.php?latex=\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}'><img src='http://latex.codecogs.com/png.latex?\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}%.png' title='\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}'><img src='http://latex.codecogs.com/png.latex?\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}%.png' title='\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}' /></a>.

The original Nose-Hoover thermostat is only valid for ergodic systems. For application to a non-ergodic system, a Nose-Hoover chain must be used ([Martyna1992](References.md)). This method introduces additional degrees of freedom that couple to the original thermostat in order to shake up the system. The length of this chain is a variable which must be determined through testing.

For path integral molecular dynamics this amounts to coupling each bead on each particle with a 1D Nose-Hoover chain. Though this may seem like a massive overhead just to achieve a canonical distribution, the time required to thermalize is usually still much less than the time required for force calculations ([Tuckerman2010](References.md)).

When integrating the Nose-Hoover equations of motion (usually through use of the Liouville operator), it is convenient to factorize the action further than the traditional Trotter expansion. To do so, we use the Suzuki-Yoshida scheme which applies a specific numerical weight to each term in the factorization where the sum of all the weights is necessarily equal to 1. This weight depends on the order of the scheme chosen, where a higher order usually results in a higher accuracy. Once an order is chosen, the weights can be determined numerically, though we simply quote the results from previous work ([Tuckerman2010](References.md)). We also determine the optimal order through testing.

Finally, one could imagine thermalization acting on a faster time scale than e.g. long-range interaction. Because of this, it is often convenient to use a smaller time step when integrating the Nose-Hoover equations of motion. This integration step is then looped over until the total time equals the time step for the rest of the simulation, <a href='http://www.codecogs.com/eqnedit.php?latex=\delta t = \Delta t / n'><img src='http://latex.codecogs.com/png.latex?\delta t = \Delta t / n%.png' title='\delta t = \Delta t / n' /></a>. Again, the optimal factor must be determined through testing.


---

**NEXT SECTION**: [Langevin\_Thermostat](Langevin_Thermostat.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)