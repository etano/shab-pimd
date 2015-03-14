**Table of Contents:**


---


# Introduction #

## Motivation ##

Often times simulations of atomic systems treat the nuclei of the atoms as point particles and make use of the Born-Oppenheimer approximation. The Born-Oppenheimer approximation exploits the fact that the nucleus of an atom is many times larger than electrons, and thus treats the nuclei as having fixed electronic degrees of freedom. A better approximation might use an ab initio technique such as Density Functional Theory (DFT) or Quantum Monte Carlo (QMC) to relax this condition. However, these approximations normally do not account for thermal expansion of a crystal lattice and can be insufficient when describing systems involving atoms with small nuclei [(Marx1995)](References.md). Conventional ab initio molecular dynamics may take into account such thermal effects, but it still treats the nuclei as classical particles, thus losing often essential quantum mechanical behavior. Thus our attention is focused on the ab initio path integral molecular dynamics method which accounts for all these effects from first principles [(Zhang2011)](References.md).
For simplicity we are working without the DFT calculation so far and use fixed potential interaction (i.e. free particles or Lennard-Jones).

## Path Integrals ##

Feynman path integrals turn a quantum mechanical problem into a series of classical problems. The _action_ is defined as

<a href='http://www.codecogs.com/eqnedit.php?latex=S=\int Ldt'><img src='http://latex.codecogs.com/png.latex?S=\int Ldt%.png' title='S=\int Ldt' /></a> ,

where _L_ is the Lagrangian. Classically, a particle follows a single path from its starting point _x_ to a new position _x'_, corresponding to the path of least action [(Peskin1995)](References.md). Using the path integral representation of quantum mechanics the particle, the particle samples all paths from _x_ to _x'_ simultaneously, where each path has a given amplitude. The probability that a particle at position _x_ will later be observed at position _x'_ is the square modulus of the sum of these interfering amplitudes.

## Path Integral Molecular Dynamics ##

The canonical one particle partition function can be written as,

<a href='http://www.codecogs.com/eqnedit.php?latex=Z(\beta)=Tr(e^{-\beta H})=\int dx \left \langle x|e^{-\beta H}|x \right \rangle'><img src='http://latex.codecogs.com/png.latex?Z(\beta)=Tr(e^{-\beta H})=\int dx \left \langle x|e^{-\beta H}|x \right \rangle%.png' title='Z(\beta)=Tr(e^{-\beta H})=\int dx \left \langle x|e^{-\beta H}|x \right \rangle' /></a>.

The Hamiltonian in the exponent is the sum of two non-commutative operators, the kinetic and potential energy. Therefore, to separate the term into the product of two non-coupled exponentials, the Trotter expansion must be used

<a href='http://www.codecogs.com/eqnedit.php?latex=e^{\lambda (A@plus;B)}=lim_{P\rightarrow \infty}[e^{\frac{\lambda}{2P}A}e^{\frac{\lambda}{P}B}e^{\frac{\lambda}{2P}A}]'><img src='http://latex.codecogs.com/png.latex?e^{\lambda (A+B)}=lim_{P\rightarrow \infty}[e^{\frac{\lambda}{2P}A}e^{\frac{\lambda}{P}B}e^{\frac{\lambda}{2P}A}]^{P}%.png' title='e^{\lambda (A+B)}=lim_{P\rightarrow \infty}[e^{\frac{\lambda}{2P}A}e^{\frac{\lambda}{P}B}e^{\frac{\lambda}{2P}A}]^{P}' /></a>.

This expression is only exact in the limit _P_ goes to infinity. The error in the Trotter expansion is order <a href='http://www.codecogs.com/eqnedit.php?latex=\tau^{2}'><img src='http://latex.codecogs.com/png.latex?\tau^{2}%.png' title='\tau^{2}' /></a> where <a href='http://www.codecogs.com/eqnedit.php?latex=\tau \equiv \beta /P'><img src='http://latex.codecogs.com/png.latex?\tau \equiv \beta /P%.png' title='tau \equiv \beta /P' /></a>. With this approximation, one is lead to the following expression for the partition function

<a href='http://www.codecogs.com/eqnedit.php?latex=Z_{P}(\beta))=(\frac{mP}{2 \pi \beta \hbar^{2}})^{\frac{P}{2}}\int dx_{1}\cdots dx_{P}exp\{-\sum_{i=1}^{P}[\frac{mP}{2 \beta \hbar^{2}}(x_{i@plus;1}-x_{i})^{2}@plus;\frac{\beta}{P}\phi(x_{i})]\}_{x_{P@plus;1}=x_{1}}'><img src='http://latex.codecogs.com/png.latex?Z_{P}(\beta))=(\frac{mP}{2 \pi \beta \hbar^{2}})^{\frac{P}{2}}\int dx_{1}\cdots dx_{P}exp\{-\sum_{i=1}^{P}[\frac{mP}{2 \beta \hbar^{2}}(x_{i+1}-x_{i})^{2}+\frac{\beta}{P}\phi(x_{i})]\}_{x_{P+1}=x_{1}}%.png' title='Z_{P}(\beta))=(\frac{mP}{2 \pi \beta \hbar^{2}})^{\frac{P}{2}}\int dx_{1}\cdots dx_{P}exp\{-\sum_{i=1}^{P}[\frac{mP}{2 \beta \hbar^{2}}(x_{i+1}-x_{i})^{2}+\frac{\beta}{P}\phi(x_{i})]\}_{x_{P+1}=x_{1}}' /></a>.

This expression is known as the _discretized path integral_. This Hamiltonian looks precisely like a set of _P_ coupled oscillators in some external potential <a href='http://www.codecogs.com/eqnedit.php?latex=U_{eff}'><img src='http://latex.codecogs.com/png.latex?U_{eff}}%.png' title='U_{eff}' /></a> where,

<a href='http://www.codecogs.com/eqnedit.php?latex=U_{eff}(x_{1},\ldots,x_{p})=\sum_{i=1}^{P}[\frac{1}{2}\frac{mP}{ \hbar^{2}}(x_{i@plus;1}-x_{i})^{2}@plus;\frac{1}{P}\phi(x_{i})]_{x_{P@plus;1}=x_{1}}'><img src='http://latex.codecogs.com/png.latex?U_{eff}(x_{1},\ldots,x_{p})=\sum_{i=1}^{P}[\frac{1}{2}\frac{mP}{ \hbar^{2} \beta^{2}}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_{i})]_{x_{P+1}=x_{1}}%.png' title='U_{eff}(x_{1},\ldots,x_{p})=\sum_{i=1}^{P}[\frac{1}{2}\frac{mP}{ \hbar^{2}}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_{i})]_{x_{P+1}=x_{1}}' /></a>

However, in order to perform a molecular dynamics simulation on this Hamiltonian, we must add in conjugate momenta. Doing so makes the partition function,

<a href='http://www.codecogs.com/eqnedit.php?latex=Z_{P}(\beta))=N \int dp_{1}\cdots dp_{P} \int dx_{1}\cdots dx_{P}exp\{-\beta [\sum_{i=1}^{P} \frac{p_{i}^{2}}{2\tilde{m}_{i}}@plus;U_{eff}(x_{1}\ldots x_{P})]\}'><img src='http://latex.codecogs.com/png.latex?Z_{P}(\beta))=N \int dp_{1}\cdots dp_{P} \int dx_{1}\cdots dx_{P}exp\{-\beta [\sum_{i=1}^{P} \frac{p_{i}^{2}}{2\tilde{m}_{i}}+U_{eff}(x_{1}\ldots x_{P})]\}%.png' title='Z_{P}(\beta))=N \int dp_{1}\cdots dp_{P} \int dx_{1}\cdots dx_{P}exp\{-\beta [\sum_{i=1}^{P} \frac{p_{i}^{2}}{2\tilde{m}_{i}}+U_{eff}(x_{1}\ldots x_{P})]\}' /></a>,

where _N_ is a constant. Notice that since these additional degrees of freedom are quadratic and decoupled (gaussian),  and thus merely change the overall normalization factor, not the overall dynamics [(Tuckerman2010)](References.md). The effective Hamiltonian becomes,

<a href='http://www.codecogs.com/eqnedit.php?latex=H(p,x)=\sum_{i=1}^{P}\frac{p_{i}^{2}}{2\bar{m}_{i}}@plus;U_{eff}(x_{1},\ldots ,x_{P})'><img src='http://latex.codecogs.com/png.latex?H(p,x)=\sum_{i=1}^{P}\frac{p_{i}^{2}}{2\bar{m}_{i}}+U_{eff}(x_{1},\ldots ,x_{P})%.png' title='H(p,x)=\sum_{i=1}^{P}\frac{p_{i}^{2}}{2\bar{m}_{i}}+U_{eff}(x_{1},\ldots ,x_{P})' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}@plus;\frac{1}{2}m\omega _{p}^{2}(x_{i@plus;1}-x_{i})^{2}@plus;\frac{1}{P}\phi(x_i)]_{x_{P@plus;1}=x_{P}}'><img src='http://latex.codecogs.com/gif.latex?H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}%.png' title='H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}' /></a>.

We are now ready to perform molecular dynamics simulations. There are now P "beads" per quantum mechanical particle. Each bead can be treated as a classical particle in a molecular dynamics simulation with the interaction Ueff. Note that every bead interacts only with the two beads that are nearest neighbors in imaginary time. Hence the one particle is represented by a necklace with P beds. To do the MD simultation, we simply implement the usual Velocity Verlet algorithm acting on each bead of every particle,

<a href='http://www.codecogs.com/eqnedit.php?latex=x(\Delta t)=x(0)@plus;\Deltat\frac{p(0)}{m}@plus;\frac{\Delta t^{2}}{2m}F(x(0))'><img src='http://latex.codecogs.com/gif.latex?x(\Delta t)=x(0)+\Delta t\frac{p(0)}{m}+\frac{\Delta t^{2}}{2m}F(x(0))%.png' title='x(\Delta t)=x(0)+\Delta t\frac{p(0)}{m}+\frac{\Delta t^{2}}{2m}F(x(0))' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=p(\Delta t)=p(0)@plus;\Delta t\frac{p(0)}{2}[F(x(0))@plus;F(x(\Delta t))]'><img src='http://latex.codecogs.com/gif.latex?p(\Delta t)=p(0)+\Delta t\frac{p(0)}{2}[F(x(0))+F(x(\Delta t))]%.png' title='p(\Delta t)=p(0)+\Delta t\frac{p(0)}{2}[F(x(0))+F(x(\Delta t))]' /></a>

In practice there are some issues that need to be addressed so that the MD simulation actually runs well. These issues will be addressed in the following sections.

# Variable Transformations #

Recall the Hamiltonian,

<a href='http://www.codecogs.com/eqnedit.php?latex=H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}@plus;\frac{1}{2}m\omega _{p}^{2}(x_{i@plus;1}-x_{i})^{2}@plus;\frac{1}{P}\phi(x_i)]_{x_{P@plus;1}=x_{P}}'><img src='http://latex.codecogs.com/gif.latex?H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}%.png' title='H(p,x)=\sum_{i=1}^{P}[\frac{p_{i}^{2}}{2\tilde{m}_i}+\frac{1}{2}m\omega _{p}^{2}(x_{i+1}-x_{i})^{2}+\frac{1}{P}\phi(x_i)]_{x_{P+1}=x_{P}}' /></a>

where

<a href='http://www.codecogs.com/eqnedit.php?latex=\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}'><img src='http://latex.codecogs.com/png.latex?\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}%.png' title='\omega_{P}=\frac{\sqrt{P}}{\beta \hbar}' /></a>.

The Trotter expansion leads to an approximation that approaches a quantum system as the number of beads goes to infinity. However, as the number of beads, _P_, is increased <a href='http://www.codecogs.com/eqnedit.php?latex=\omega_{P}'><img src='http://latex.codecogs.com/png.latex?\omega_{P}^{2}%.png' title='\omega_{P}^{2}' /></a> scales linearly with the number of beads. Thus the entire second term scales linearly with the number of beads. Comparatively the first term does not depend on _P_ and the external potential term varies inversely with _P_. As a result, the forces due to the harmonic term, which couples all of the beads, dominate. This strong coupling restricts the motion of the necklace as well as the size of the time step that may be used. These challenges are addressed with a change of variables and a time scale iteration algorithm that incorporates two different time scales [(Tuckerman2010)](References.md).

In order to avoid these issues, we wish to uncouple the harmonic term through a variable transformation. In this project we investigate the effects of two different approaches: the Staging transformation and the normal mode transformation.

## Staging Transformation ##

The first change of variables discussed is known as a _staging transformation_,

<a href='http://www.codecogs.com/eqnedit.php?latex=x_{1}=u_{1}'><img src='http://latex.codecogs.com/gif.latex?x_{1}=u_{1}%.png' title='x_{1}=u_{1}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=x_{i}=u_{i}@plus;\frac{i-1}{i}x_{i@plus;1}@plus;\frac{1}{i}u_{1}'><img src='http://latex.codecogs.com/gif.latex?x_{i}=u_{i}+\frac{i-1}{i}x_{i+1}+\frac{1}{i}u_{1}%.png' title='x_{i}=u_{i}+\frac{i-1}{i}x_{i+1}+\frac{1}{i}u_{1}' /></a> , where <a href='http://www.codecogs.com/eqnedit.php?latex=i=2,\ldots, P'><img src='http://latex.codecogs.com/gif.latex?i=2,\ldots, P%.png' title='i=2,\ldots, P' /></a>.

With this transformation the masses are no longer uniform. Instead they become,

<a href='http://www.codecogs.com/eqnedit.php?latex=m_{1}=0'><img src='http://latex.codecogs.com/png.latex?m_{1}=0%.png' title='m_{1}=0' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=m_{k}=\frac{k}{k-1}m,'><img src='http://latex.codecogs.com/png.latex?m_{k}=\frac{k}{k-1}m,%.png' title='m_{k}=\frac{k}{k-1}m,' /></a>   where    <a href='http://www.codecogs.com/eqnedit.php?latex=k=2,\ldots,P'><img src='http://latex.codecogs.com/png.latex?k=2,\ldots,P.%.png' title='k=2,\ldots,P' /></a>

so that the new masses become <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{m}_{1}=m'><img src='http://latex.codecogs.com/png.latex?\tilde{m}_{1}=m%.png' title='\tilde{m}_{1}=m' /></a> and <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{m}_{k}=\frac{k}{k-1}m'><img src='http://latex.codecogs.com/png.latex?\tilde{m}_{k}=\frac{k}{k-1}m%.png' title='\tilde{m}_{k}=\frac{k}{k-1}m' /></a> .

Furthermore, this transformation leads to the following equations of motion

<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}%.png' title='\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}%.png' title='\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}' /></a>

Finally using the chain rule we arrive upon the following forces

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\frac{1}{P}\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}@plus;\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}+\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}+\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]' /></a>.

## Normal Mode Transformation ##

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

# Thermostats #

As written in the introduction, our molecular dynamics simulation is in the constant energy micro-canonical ensemble. To sample the partition function at a constant temperature canonical ensemble a thermostat is usually required. Note that necessarily this is unphysical since one is adding degrees of freedom to the system which are not truly present and changing the conserved quantity (the Hamiltonian). Nevertheless, there has been much success with correctly reproducing the correct canonical distribution through thermostatting. We compare two methods: The Nose-Hoover thermostat deterministically introduces new degrees of freedom that lead to a constant temperature. It is considered the _gold standard_ of thermostatting. The Langevin thermostat makes use of the Langevin Equation to mimic a heat bath and it thus stochastic. It employs the Einstein relation between friction and the random force to set the temperature.

## Nose-Hoover Thermostat ##

The Nose-Hoover thermostat centers on an extended system method ([Hunenberger2005](References.md)). The basic Nose-Hoover thermostat introduces one additional degree of freedom, <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{s}'><img src='http://latex.codecogs.com/png.latex?\tilde{s}%.png' title='\tilde{s}' /></a>, such that the Lagrangian becomes

<a href='http://www.codecogs.com/eqnedit.php?latex=L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})@plus;\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})'><img src='http://latex.codecogs.com/png.latex?L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})+\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})%.png' title='L_{Nose}(\tilde{r},\dot{\tilde{r}},\tilde{s},\dot{\tilde{s}})=\frac{1}{2}\sum_{i=1}^{N}m_{i}\tilde{s}^{2}\dot{\tilde{r}}_{i}^{2}-U(\tilde{r})+\frac{1}{2}Q\dot{\tilde{s}}^{2}-gk_{B}T_{0}ln(\tilde{s})' /></a>.

where the ~ denotes extended system variables. Now our conserved quantity is simply the original Hamiltonian plus the energy of the added Nose-Hoover degrees of freedom. It follows that the equations of motion are

<a href='http://www.codecogs.com/eqnedit.php?latex=\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}'><img src='http://latex.codecogs.com/png.latex?\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}%.png' title='\ddot{\tilde{r}}_{i}=\frac{F_{i}}{m_{i}\tilde{s}^2}-\frac{2\dot{\tilde{s}}\dot{\tilde{r}}_{i}}{\tilde{s}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}'><img src='http://latex.codecogs.com/png.latex?\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}%.png' title='\ddot{\tilde{s}}_{i}=\frac{\sum_{i=1}^{N}m_{i}\dot{\tilde{r}}_{i}^2\tilde{s}^2-gk_{B}T_{0}}{Q\tilde{s}}' /></a>.

The original Nose-Hoover thermostat is only valid for ergodic systems. For application to a non-ergodic system, a Nose-Hoover chain must be used ([Martyna1992](References.md)). This method introduces additional degrees of freedom that couple to the original thermostat in order to shake up the system. The length of this chain is a variable which must be determined through testing.

For path integral molecular dynamics this amounts to coupling each bead on each particle with a 1D Nose-Hoover chain. Though this may seem like a massive overhead just to achieve a canonical distribution, the time required to thermalize is usually still much less than the time required for force calculations ([Tuckerman2010](References.md)).

When integrating the Nose-Hoover equations of motion (usually through use of the Liouville operator), it is convenient to factorize the action further than the traditional Trotter expansion. To do so, we use the Suzuki-Yoshida scheme which applies a specific numerical weight to each term in the factorization where the sum of all the weights is necessarily equal to 1. This weight depends on the order of the scheme chosen, where a higher order usually results in a higher accuracy. Once an order is chosen, the weights can be determined numerically, though we simply quote the results from previous work ([Tuckerman2010](References.md)). We also determine the optimal order through testing.

Finally, one could imagine thermalization acting on a faster time scale than e.g. long-range interaction. Because of this, it is often convenient to use a smaller time step when integrating the Nose-Hoover equations of motion. This integration step is then looped over until the total time equals the time step for the rest of the simulation, <a href='http://www.codecogs.com/eqnedit.php?latex=\delta t = \Delta t / n'><img src='http://latex.codecogs.com/png.latex?\delta t = \Delta t / n%.png' title='\delta t = \Delta t / n' /></a>. Again, the optimal factor must be determined through testing.

## Langevin Thermostat ##

An alternative to the deterministic Nose Hoover thermostat is the stochastic Langevin thermostat ([Ceriotti2010](References.md)). This thermostat modifies the forces by adding a friction term and a thermal noise term

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}'><img src='http://latex.codecogs.com/png.latex?\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}%.png' title='\frac{\partial}{\partial t}\tilde{q}_{i}^{(k)}=\frac{\tilde{p}_{i}^{(k)}}{m_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}@plus;\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)'><img src='http://latex.codecogs.com/png.latex?\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}+\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)%.png' title='\frac{\partial}{\partial t}\tilde{p}_{i}^{(k)}=-m_{i}\omega_k^{2}\tilde{q}_{i}^{(k)}-\gamma^{(k)}\tilde{p}_{i}^{(k)}+\sqrt{\frac{2m_{i}\gamma^{(k)}}{\beta_{n}}}\xi_{i}^{(k)}(t)' /></a> ,

where <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> are the friction coefficients and <a href='http://www.codecogs.com/eqnedit.php?latex=\xi_{i}^{(k)}(t)'><img src='http://latex.codecogs.com/png.latex?\xi_{i}^{(k)}(t)%.png' title='\xi_{i}^{(k)}(t)' /></a> is a normally-distributed random force. Here our conserved quantity changes from the total energy of the system to the total energy minus the the accumulated heat absorbed by the thermostat. The Langevin thermostat ensures canonical dynamics no matter the value of <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a>, though better values should lead to faster convergence. The optimal values of <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> can be determined recursively ([Ceriotti2010](References.md)). However, for our purposes, this made little difference to the efficiency of the simulation. Thus we let each <a href='http://www.codecogs.com/eqnedit.php?latex=\gamma^{(k)}'><img src='http://latex.codecogs.com/png.latex?\gamma^{(k)}%.png' title='\gamma^{(k)}' /></a> be unity.

# Results #

Below we discuss all our results for both non-interacting and interacting systems. Note that for all simulations we are in 1D and contained in a harmonic well external potential.

## Nose-Hoover Testing ##

As mentioned previously, in order to fairly compare the Nose-Hoover thermostat, it was necessary to find the optimal parameters for the length of the Nose-Hoover chain, the order of the Suzuki-Yoshida expansion, and the number of discretized time steps required for proper thermalization. In order to do so, we compare convergence efficiencies for a single particle in a harmonic trap for a series of various parameters. Specifically, we look at a low temperature, <a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> =10.0, and a high temperature, <a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0, in order to find settings which are suitable for a range of temperatures.

<a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0<img src='http://www.etano.net/phys466/EffBarVEvBeta-NHT10.png' />

<a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0
<img src='http://www.etano.net/phys466/EffBarVEvBeta-NHT1.png' />

As one can clearly see, at low temperatures, the settings (2,4,1) are most efficient. This corresponds to a Nose-Hoover length of 2, a Suzuki-Yoshida expansion of order 4, and a single thermalization step. However, at high temperatures these settings result in zero efficiency. This implies that a single thermalization step is not adequate for proper coupling to the heat bath. Thus, we settle upon the second set of settings (2,4,2) for the remainder of our methods comparison. Note that for all comparisons we choose the virial energy estimator which proved more reliable than the primitive estimator.

## Non-Interacting System ##

To study the non-interacting system we examine only a single particle in one dimension in a harmonic trap. (Studying more particles without interaction would prove unnecessary.) Below we plot convergence in the number of time slices for each method combination. Note that we are showing the scaled virial energy estimator. It is being scaled by the analytical classical energy, making 0 the classical energy. The black line represents the exact quantum mechanical energy found through exact diagonalization. For all simulations, we used 1000000 sweeps with a time step of 0.01 and a block size of 10000.

**No Transformation, Nose-Hoover:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_010.png' />

**No Transformation, Langevin:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_020.png' />

**Staging, Nose-Hoover:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_110.png' />

**Staging, Langevin:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_120.png' />

**Normal Mode, Nose-Hoover:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_210.png' />

**Normal Mode, Langevin:**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_220.png' />

It is difficult to tell how each method performs relative to the others from these plots. However, it does appear that methods using the Langevin thermostat are better converged as well as are those using some variable transformation. For all the methods we notice an anomalous gap between the converged result and the exact quantum mechanical result around <a href='http://www.codecogs.com/eqnedit.php?latex=\beta = 4'><img src='http://latex.codecogs.com/png.latex?\beta = 4%.png' title='\beta = 4' /></a>. It has yet to be confirmed, but this feature is likely a cause of the non-ergodicity of the harmonic trap. As mentioned in a previous section, any sort of thermostatting is in essence _unphysical_ since it adds degrees of freedom not necessarily found in the actual system. What is important, however is that all the methods are converging to approximately the same result making their direct comparison a viable option.

### Comparison of Methods ###

For a better comparison of the methods we must look at the relative simulation durations and errors. Using these two numbers, we can then compute the efficiency of each simulation, <a href='http://www.codecogs.com/eqnedit.php?latex=\zeta \sim 1/(T \sigma^{2})'><img src='http://latex.codecogs.com/png.latex?\zeta \sim 1/(T \sigma^{2})%.png' title='\zeta \sim 1/(T \sigma^{2})' /></a>.

**Consistency Between Methods**
<img src='http://www.etano.net/phys466/VEScaledvBeta-beads_m32_Eff.png' />

This plot demonstrates the each method combination converges to approximately the same result. Here the number of times slices is 32. Thus, a comparison of the methods must look at second order results such as the error and subsequently, the efficiency.

**Error of Methods**
<img src='http://www.etano.net/phys466/ErrVEvBeta-beads_m32_Eff.png' />

Errors between the methods are fairly similar for most temperatures. However, at high temperatures the simulations without a variable transformation perform worse. This is because, as mentioned previously, without a transformation all modes are coupled. At high temperatures this causes the particles to oscillate more wildly.

**Time per Block of Methods**
<img src='http://www.etano.net/phys466/tvBeta-beads_m32_Eff.png' />

The above plot compares the average time required to complete a single block for each method combination. The two quickest methods are using the Langevin thermostat, one without a transformation and one with the staging transformation. Methods using the normal mode transformation tend to be slower. This is most likely because the orthogonal matrix multiplication required for the method must loop through all time slices. On the other hand, the staging transformation is defined by a recursion relation requiring information on the previous time slice. This discrepancy could be made smaller using Fast Fourier Transforms with is order N log(N). However, this is still larger than the order of the staging transformation.

**Efficiency of Methods**
<img src='http://www.etano.net/phys466/EffVEvBeta-beads_m32_Eff.png' />

As one can see, the combination of a staging transformation coupled to a Langevin thermostat proved to be the most efficient method for swift convergence. The discrepancies between the other methods are debatable, though it seems that having no transformation at all is never favorable.

## Interacting System ##

Now we turned to the more interesting case of interacting quantum particle. We chose particles trapped in a harmonic potential with a Lennard Jones interaction (Epsilon = 1, [r0](https://code.google.com/p/shab-pimd/source/detail?r=0) = 1) . To really do ab initio calculations one would need to implement a running DFT calculation at each time slice. This could be added to the code but we worked with a fixed potential for simplicity.

### Two-Particle Check ###

To have an idea if the results we are getting from the simulation actually approach the correct quantum mechanical result, we started with only two particles. The plot shows again the total energy (by the Virial Energy Estimator VE) at different temperatures beta and with different numbers of beads per particle. For comparison we added the classical result from statistical mechanics (green) obtained by numerical integration of the partition function and the quantum result (blue) obtained by direct diagonalization of a discretized Hamiltonian. As in the non-interacting case the lines with different number beads of per particle interpolate between these two limits. But convergence to the quantum behavior is much slower. Even for 2 times 64 beads, there is still a noticeable offset from the quantum result. For only one bead per particle, which should be a result classical Molecular Dynamics simulation we do not understand yet the reason for the offset from the result by integration of the partition function.

<img src='http://www.etano.net/phys466/VEvBeta-beads2_121b.png' />
<img src='http://www.etano.net/phys466/BSvBeta-beads2_121b.png' />

### Five-Particle Study ###
After we have been convinced that our code converges towards the exact quantum limit for large numbers of beads, we did again the comparison of the different transformations and thermostats. To have little more interesting system, we chose 5 particles with Lennard-Jones interaction in a harmonic trap.

#### Comparison of Methods ####

The first plot shows that all methods converge to the same result. This is for 32 beads per particle (thus 160 beads total), where we are still quite far away from the quantum limit. Still it should be possible to analyze the efficiency and draw some results from it that are also true for larger numbers of beads.

<img src='http://www.etano.net/phys466/VEvBeta-beads_int_m32_eff.png' />

So we look at the CPU-time it takes to run the different methods (with the same time step, same number of steps) and at the statistical error (uncertainty). One can see again that the Nose-Hoover thermostat really slows the calculation because of all the additional degrees of freedom. This effect should reduce however for increasing number of interacting particles because the force/energy calculation is the only part that scales as N<sup>2</sup> and will dominate for large N, so that the extra time of the [Nose\_Hoover\_Thermostat](Nose_Hoover_Thermostat.md) thermostat will be less relevant for the efficiency. So far we could not find an advantage of the Nose-Hoover over the Langevin thermostat by comparing the statistical error. As expected one can also see that the statistical error of the no-transformation calculations is a lot higher. Thus, a transformation is important to achieve fast convergence.


Finally we calculated the efficiency. For the case of 5 particles the [Langevin\_Thermostat](Langevin_Thermostat.md) thermostat with the [Staging\_Transformation](Staging_Transformation.md) clearly outperforms all other combinations. So this result agrees with the non-interacting case.

<img src='http://www.etano.net/phys466/ErrVEvBeta-beads_int_m32_eff.png' />
<img src='http://www.etano.net/phys466/tvBeta-beads_int_m32_eff.png' />


<img src='http://www.etano.net/phys466/EffVEvBeta-beads_int_m32_eff.png' />

#### Supplementary Observables ####

We also looked at some other observables like the density and the pair correlation function.

The density shows how the five particle line up in the one dimensional trap for low temperature. One in the middle and two on each side (Plot only show half the trap for x > 0 ). One can see how the system thermally expands with increasing temperature due to the non-harmonic parts of the Lennard-Jones interaction. The density was calculated via a histogram of the bead centroids.

<img src='http://www.etano.net/phys466/RDen-beads_021_m32_grr.png' />

The pair correlation (also calculated from the bead centroids) shows the
same behavior of localized particles at low temperatures, which thermally expand for higher temperatures. The decrease of the height of the peaks is due to the finite number of particles (without periodic boundaries).
<img src='http://www.etano.net/phys466/Grr-beads_021_m32_grr.png' />

# Conclusions #

We are able to use PIMD to model a quantum mechanical system. At high temperatures, PIMD simulations approach a classical result where all of the beads move close together. At lower temperatures, the beads spread out, corresponding to a quantum mechanical regime. Our simulations clearly showed the relationship between bead spreading and classical or quantum mechanical behavior.

We found that the number of time slices, or beads on the necklace, required for quantum convergence varies depending on the system. Also, the time step for high temperatures needs to be smaller than what we used to run our simulations. However, the focus of our project was on the comparison of different method combinations, so we used the same time step in every simulation.

The best combination of methods is clearly the Langevin thermostat with the staging transformation. This conclusion is based upon the time that it took to run the simulation and the errors in the values predicted by the simulation. The next best combination is the Langevin thermostat with the normal mode transformation. Even with our optimizations, the Nose-Hoover thermostat in combination with any of the transformations did not compare to the use of the Langevin thermostat with the staging transformation.

# Future Directions #

In the future, it would be interesting to apply the PIMD approach to an actual system using our optimal settings. Analysis of systems such as the interaction of hydrogen with a metal lattice would benefit from the use of the PIMD method, because it would incorporate quantum mechanical effects due to the small mass of the hydrogen atom. For this application, PIMD would be combined with an ab initio method such as DFT. For further comparisons, PIMD should be compared with other methods such as the path integral Monte Carlo method, looking at efficiency.

# References #

  1. Ceriotti, M., M. Parrinello, T.E. Markland, and D.E. Manolopoulos, "Efficient Stochastic Thermostatting of Path Integral Molecular Dynamics." Journal of Chemical Physics **133** 124104-124104-13 (2010).
  1. Hunenberger, P.H., _Thermostat Algorithms for Molecular Dynamics Simulations_. Springer Berlin / Heidelberg (2005).
  1. Martyna, G.J., M.L. Klein, and M. Tuckerman, "Nose-Hoover Chains: The Canonical Ensemble Via Continuous Dynamics." Journal of Chemical Physics **97**, 2635-2643 (1992).
  1. Marx, D. and M. Parrinello, "Ab Initio Path Integral Molecular Dynamics: Basic Ideas." Journal of Chemical Physics **104**, 4077-4082 (1995).
  1. Peskin, M.E. and D.V. Schroeder, _An Introduction of Quantum Field Theory_. Westview Press (1995).
  1. Tuckerman, M., "Path Integration via Molecular Dynamics." _Quantum Simulations of Complex Many-Body Systems: From Theory to Algorithms_, Lecture Notes **10**, 269-298 (2002).
  1. Tuckerman, M., _Statistical Mechanics: Theory and Molecular Simulation_. Oxford University Press (2010).
  1. Zhang, C. and A. Michaelides, "Quantum Nuclear Effects on the Location of Hydrogen Above and Below the Palladium (100) Surface." Surface Science **605**, 689-694 (2011).