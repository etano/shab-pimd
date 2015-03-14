Another change of variables that may be used to uncouple the harmonic term is the normal mode transformation [(Tuckerman2010)](References.md). This method may either be implemented with a decomposition into Fourier modes. Such a change corresponds to

<a href='http://www.codecogs.com/eqnedit.php?latex=x_{k}=\sqrt{P}\sum_{l=1}^{P}O_{kl}^{T}u_{l}'><img src='http://latex.codecogs.com/gif.latex?x_{k}=\sqrt{P}\sum_{l=1}^{P}O_{kl}^{T}u_{l}%.png' title='x_{k}=\sqrt{P}\sum_{l=1}^{P}O_{kl}^{T}u_{l}' /></a>,

where <a href='http://www.codecogs.com/eqnedit.php?latex=O_{kl}^{T}'><img src='http://latex.codecogs.com/gif.latex?O_{kl}^{T}%.png' title='O_{kl}^{T}' /></a> is an orthogonal matrix that diagonalizes a matrix of the form

<a href='http://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix} 2 & -1 & 0 & 0 & 0 & 0 & 0 & -1\\ -1 & 2 & -1 & 0 & 0 & 0 & 0 & 0\\ 0 & -1 & 2 & -1 & 0 & 0 & 0& 0\\ 0 & 0 & \ddots & \ddots & \ddots & & \vdots &\vdots\\ \vdots & \vdots & & \ddots & \ddots & \ddots & 0 & 0\\ 0 & 0 & 0 & 0 & -1 & 2 & -1 & 0\\ 0 & 0 & 0 & 0 & 0 & -1 & 2 & -1\\ -1 & 0 & 0 & 0 & 0 & 0 & -1 & 2 \end{bmatrix}'><img src='http://latex.codecogs.com/gif.latex?\begin{bmatrix} 2 & -1 & 0 & 0 & 0 & 0 & 0 & -1\\ -1 & 2 & -1 & 0 & 0 & 0 & 0 & 0\\ 0 & -1 & 2 & -1 & 0 & 0 & 0& 0\\ 0 & 0 & \ddots & \ddots & \ddots & & \vdots &\vdots\\ \vdots & \vdots & & \ddots & \ddots & \ddots & 0 & 0\\ 0 & 0 & 0 & 0 & -1 & 2 & -1 & 0\\ 0 & 0 & 0 & 0 & 0 & -1 & 2 & -1\\ -1 & 0 & 0 & 0 & 0 & 0 & -1 & 2 \end{bmatrix}%.png' title='\begin{bmatrix} 2 & -1 & 0 & 0 & 0 & 0 & 0 & -1\\ -1 & 2 & -1 & 0 & 0 & 0 & 0 & 0\\ 0 & -1 & 2 & -1 & 0 & 0 & 0& 0\\ 0 & 0 & \ddots & \ddots & \ddots & & \vdots &\vdots\\ \vdots & \vdots & & \ddots & \ddots & \ddots & 0 & 0\\ 0 & 0 & 0 & 0 & -1 & 2 & -1 & 0\\ 0 & 0 & 0 & 0 & 0 & -1 & 2 & -1\\ -1 & 0 & 0 & 0 & 0 & 0 & -1 & 2 \end{bmatrix}' /></a>.

For this transformation, the equations of motion are the same as before


<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}%.png' title='\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{p}_{i}=-\frac{\partial U_{eff}^{normal}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{p}_{i}=-\frac{\partial U_{eff}^{normal}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}%.png' title='\dot{p}_{i}=-\frac{\partial U_{eff}^{normal}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}' /></a>.

However, now the masses are given by

<a href='http://www.codecogs.com/eqnedit.php?latex=m_{k}=m\lambda_{k}'><img src='http://latex.codecogs.com/png.latex?m_{k}=m\lambda_{k}%.png' title='m_{k}=m\lambda_{k}' /></a>

where <a href='http://www.codecogs.com/eqnedit.php?latex=\lambda_{k}'><img src='http://latex.codecogs.com/png.latex?\lambda_{k}%.png' title='\lambda_{k}' /></a> are the eigenvalues of <a href='http://www.codecogs.com/eqnedit.php?latex=O_{kl}^{T}'><img src='http://latex.codecogs.com/gif.latex?O_{kl}^{T}%.png' title='O_{kl}^{T}' /></a>.For this variable transformation, the forces are given by

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}}'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}}%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{1}{P}\frac{\partial \phi}{\partial u_{k}}=\frac{1}{\sqrt{P}}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}O_{kl}^{T}}'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{k}}=\frac{1}{\sqrt{P}}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}O_{kl}^{T}}%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{k}}=\frac{1}{\sqrt{P}}\sum_{l=1}^{P}\frac{\partial \phi}{\partial x_{l}O_{kl}^{T}}' /></a>.



---

**NEXT SECTION**: [Thermostats](Thermostats.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)