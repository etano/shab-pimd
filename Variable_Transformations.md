# Variable Transformations #

Recall the Hamiltonian,

<a href='http://www.codecogs.com/eqnedit.php?latex=H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}@plus;\frac{1}{2}m\omega _{p}^{2}(x_{i@plus;1}-x_{i})^{2}@plus;\frac{1}{P}\phi(x_i)]_{x_{P@plus;1}=x_{P}}'><img src='http://latex.codecogs.com/gif.latex?H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}%.png' title='H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}' /></a>

where

<a href='http://www.codecogs.com/eqnedit.php?latex=\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}'><img src='http://latex.codecogs.com/png.latex?\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}%.png' title='\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}' /></a>.

The Trotter expansion leads to an approximation that approaches a quantum system as the number of beads goes to infinity. However, as the number of beads, _P_, is increased <a href='http://www.codecogs.com/eqnedit.php?latex=\omega_{P}'><img src='http://latex.codecogs.com/png.latex?\omega_{P}^{2}%.png' title='\omega_{P}^{2}' /></a> scales linearly with the number of beads. Thus the entire second term scales linearly with the number of beads. Comparatively the first term does not depend on _P_ and the external potential term varies inversely with _P_. As a result, the forces due to the harmonic term, which couples all of the beads, dominate. This strong coupling restricts the motion of the necklace as well as the size of the time step that may be used. These challenges are addressed with a change of variables and a time scale iteration algorithm that incorporates two different time scales [(Tuckerman2010)](References.md).

In order to avoid these issues, we wish to uncouple the harmonic term through a variable transformation. In this project we investigate the effects of two different approaches:
  * [Staging\_Transformation](Staging_Transformation.md)
  * [Normal\_Mode\_Transformation](Normal_Mode_Transformation.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)