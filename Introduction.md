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

In practice there are some issues that need to be addressed so that the MD simulation actually runs well. These issues will be addressed in the following sections [Variable\_Transformations](Variable_Transformations.md) and [Thermostats](Thermostats.md).


---

**NEXT SECTION**: [Variable\_Transformations](Variable_Transformations.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)