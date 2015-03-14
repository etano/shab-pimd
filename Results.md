Below we discuss all our results for both non-interacting and interacting systems. Note that for all simulations we are in 1D and contained in a harmonic well external potential.

---

**Table of Contents:**


---

# Nose-Hoover Testing #

As mentioned previously, in order to fairly compare the Nose-Hoover thermostat, it was necessary to find the optimal parameters for the length of the Nose-Hoover chain, the order of the Suzuki-Yoshida expansion, and the number of discretized time steps required for proper thermalization. In order to do so, we compare convergence efficiencies for a single particle in a harmonic trap for a series of various parameters. Specifically, we look at a low temperature, <a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> =10.0, and a high temperature, <a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0, in order to find settings which are suitable for a range of temperatures.

<a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0<img src='http://www.etano.net/phys466/EffBarVEvBeta-NHT10.png' />

<a href='http://www.codecogs.com/eqnedit.php?latex=\beta'><img src='http://latex.codecogs.com/png.latex?\beta%.png' title='\beta' /></a> = 1.0
<img src='http://www.etano.net/phys466/EffBarVEvBeta-NHT1.png' />

As one can clearly see, at low temperatures, the settings (2,4,1) are most efficient. This corresponds to a Nose-Hoover length of 2, a Suzuki-Yoshida expansion of order 4, and a single thermalization step. However, at high temperatures these settings result in zero efficiency. This implies that a single thermalization step is not adequate for proper coupling to the heat bath. Thus, we settle upon the second set of settings (2,4,2) for the remainder of our methods comparison. Note that for all comparisons we choose the virial energy estimator which proved more reliable than the primitive estimator.

# Non-Interacting Results #

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

## Comparison of Methods ##

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

# Results for Interacting Particles #

Now we turned to the more interesting case of interacting quantum particle. We chose particles trapped in a harmonic potential with a Lennard Jones interaction (Epsilon = 1, [r0](https://code.google.com/p/shab-pimd/source/detail?r=0) = 1) . To really do ab initio calculations one would need to implement a running DFT calculation at each time slice. This could be added to the code but we worked with a fixed potential for simplicity.

## Two-Particle Check ##

To have an idea if the results we are getting from the simulation actually approach the correct quantum mechanical result, we started with only two particles. The plot shows again the total energy (by the Virial Energy Estimator VE) at different temperatures beta and with different numbers of beads per particle. For comparison we added the classical result from statistical mechanics (green) obtained by numerical integration of the partition function and the quantum result (blue) obtained by direct diagonalization of a discretized Hamiltonian. As in the non-interacting case the lines with different number beads of per particle interpolate between these two limits. But convergence to the quantum behavior is much slower. Even for 2 times 64 beads, there is still a noticeable offset from the quantum result. For only one bead per particle, which should be a result classical Molecular Dynamics simulation we do not understand yet the reason for the offset from the result by integration of the partition function.

<img src='http://www.etano.net/phys466/VEvBeta-beads2_121b.png' />
<img src='http://www.etano.net/phys466/BSvBeta-beads2_121b.png' />

## Five-Particle Study ##
After we have been convinced that our code converges towards the exact quantum limit for large numbers of beads, we did again the comparison of the different transformations and thermostats. To have little more interesting system, we chose 5 particles with Lennard-Jones interaction in a harmonic trap.

### Comparison of Methods ###

The first plot shows that all methods converge to the same result. This is for 32 beads per particle (thus 160 beads total), where we are still quite far away from the quantum limit. Still it should be possible to analyze the efficiency and draw some results from it that are also true for larger numbers of beads.

<img src='http://www.etano.net/phys466/VEvBeta-beads_int_m32_eff.png' />

So we look at the CPU-time it takes to run the different methods (with the same time step, same number of steps) and at the statistical error (uncertainty). One can see again that the Nose-Hoover thermostat really slows the calculation because of all the additional degrees of freedom. This effect should reduce however for increasing number of interacting particles because the force/energy calculation is the only part that scales as N<sup>2</sup> and will dominate for large N, so that the extra time of the [Nose\_Hoover\_Thermostat](Nose_Hoover_Thermostat.md) thermostat will be less relevant for the efficiency. So far we could not find an advantage of the Nose-Hoover over the Langevin thermostat by comparing the statistical error. As expected one can also see that the statistical error of the no-transformation calculations is a lot higher. Thus, a transformation is important to achieve fast convergence.


Finally we calculated the efficiency. For the case of 5 particles the [Langevin\_Thermostat](Langevin_Thermostat.md) thermostat with the [Staging\_Transformation](Staging_Transformation.md) clearly outperforms all other combinations. So this result agrees with the non-interacting case.

<img src='http://www.etano.net/phys466/ErrVEvBeta-beads_int_m32_eff.png' />
<img src='http://www.etano.net/phys466/tvBeta-beads_int_m32_eff.png' />


<img src='http://www.etano.net/phys466/EffVEvBeta-beads_int_m32_eff.png' />

### Supplementary Observables ###

We also looked at some other observables like the density and the pair correlation function.

The density shows how the five particle line up in the one dimensional trap for low temperature. One in the middle and two on each side (Plot only show half the trap for x > 0 ). One can see how the system thermally expands with increasing temperature due to the non-harmonic parts of the Lennard-Jones interaction. The density was calculated via a histogram of the bead centroids.

<img src='http://www.etano.net/phys466/RDen-beads_021_m32_grr.png' />

The pair correlation (also calculated from the bead centroids) shows the
same behavior of localized particles at low temperatures, which thermally expand for higher temperatures. The decrease of the height of the peaks is due to the finite number of particles (without periodic boundaries).
<img src='http://www.etano.net/phys466/Grr-beads_021_m32_grr.png' />


---

**NEXT SECTION**: [Conclusions](Conclusions.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)