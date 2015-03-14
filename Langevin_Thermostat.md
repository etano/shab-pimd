# Langevin Thermostat #

An alternative to the deterministic Nose Hoover thermostat is the stochastic Langevin thermostat ([Ceriotti2010](References.md)). This thermostat modifies the forces by adding a friction term and a thermal noise term

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}'><img src='http://latex.codecogs.com/png.latex?\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}%.png' title='\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}@plus;\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)'><img src='http://latex.codecogs.com/png.latex?\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}+\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)%.png' title='\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}+\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)' /></a> ,

where <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> are the friction coefficients and <a href='http://www.codecogs.com/eqnedit.php?latex=\xi_{i}^{(k)}(t)'><img src='http://latex.codecogs.com/png.latex?\xi_{i}^{(k)}(t)%.png' title='\xi_{i}^{(k)}(t)' /></a> is a normally-distributed random force. Here our conserved quantity changes from the total energy of the system to the total energy minus the the accumulated heat absorbed by the thermostat. The Langevin thermostat ensures canonical dynamics no matter the value of <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a>, though better values should lead to faster convergence. The optimal values of <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> can be determined recursively ([Ceriotti2010](References.md)). However, for our purposes, this made little difference to the efficiency of the simulation. Thus we let each <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> be unity.


---

**NEXT SECTION**: [Results](Results.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)