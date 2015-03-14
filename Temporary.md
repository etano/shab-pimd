# Results for Interacting Particles #

Now we turned to the more interisting case of interacting quantum particle. We chose particles trapped in a harmonic potential with a Lennard Jones interaction (Epsilon = 1, [r0](https://code.google.com/p/shab-pimd/source/detail?r=0) = 1) . To really do ab initio calculations one would need to implement a running DFT calculation at each time slice. This could be added to the code but we worked with a fixed potential for simplicity.

## Two-Particle Check ##

To have an idea if the results we are getting from the similation actually approach the correct quantum mechanical result, we started with only two particles. The plot shows again the total energy (by the Virial Energy Estimator VE) at different temperatures beta and with different numbers of beads per particle. For comparisson we added the classical result from statistical mechanics (green) obtained by numberical integration of the partition function and the quantum result (blue) obtained by direct dioganilzation of a discretized Hamiltonian. As in the non-interacting case the lines with different number beads of per particle interpolate between these two limits. But convergence to the quantum behaviour is much slower. Even for 2 times 64 beads, there is still a noticable offset from the quantum result. For onyl one bead per particle, which should be a regualt classical Molecular Dynamics simulation we do not understand yet the reason for the offset from the result by integration of the partition function.

<img src='http://www.etano.net/phys466/VEvBeta-beads2_121b.png' />
<img src='http://www.etano.net/phys466/BSvBeta-beads2_121b.png' />

## Five-Particle Study ##
After we have been convinced that our code convergeces towards the exact quantum limit for large numbers of beads, we did again the comparison of the different transformations and thermostats. To have little more intersting system, we chose 5 particles with Lennard-Jones interaction in a harmonic trap.

### Comparison of Methods ###

The first plot shows that all methods converge to the same result. This is for 32 beads per particle (thus 160 beads total), where we are still quite far away from the quantum limit. Still it should be possible to analyse the efficiency and draw some results from it that are also true for larger numbers of beads.

<img src='http://www.etano.net/phys466/VEvBeta-beads_int_m32_eff.png' />

So we look at the CPU-time it takes to run the different methods (with the same timestep, same number of steps) and at the statistical error (uncertianty). One can see again that the Nose-Hoover thermostat really slows the calculation because of all the additional degrees of freedom. This effect should reduce however for increasing number of interacting particles because the force/energy calculation is the only part that scales as N<sup>2</sup> and will dominate for large N, so that the extra time of the [Nose\_Hoover\_Thermostat](Nose_Hoover_Thermostat.md) thermostat will be less relevant for the efficiency. So far we could not find an advantage of the Nose-Hoover over the Langevin thermostat by comparing the statistical error. As expected one can also see that the statistical error of the no-transformation calculations is a lot higher. Thus, a transformation is important to achieve fast convergence.

<img src='http://www.etano.net/phys466/ErrVEvBeta-beads_int_m32_eff.png' />
<img src='http://www.etano.net/phys466/tvBeta-beads_int_m32_eff.png' />

Finally we calculated the efficiency. For the case of 5 particles the [Langevin\_Thermostat](Langevin_Thermostat.md) thermostat with the staging transformation clearly outperforms all other combinations. So this result agrees with the non-interacting case.

<img src='http://www.etano.net/phys466/EffVEvBeta-beads_int_m32_eff.png' />

### Supplementary Observables ###

We also looked at some other observables like the density and the pair correlation function.

The density shows how the five particle line up in the one dimensional trap for low temperature. One in the middle and two on each side (Plot only show half the trap for x > 0 ). One can see how the system thermally expands with increasing temperature due to the non-harmonic parts of the Lennard-Jones interaction. The density was calculated via a histogram of the bead centroids.

<img src='http://www.etano.net/phys466/RDen-beads_021_m32_grr.png' />

The pair corrlation (also calculated from the bead centroids) shows the
same behaviour of loalized particles at low temperatures, which thermally expand for higer temperatures. The decrease of the highth of the peaks is due to the finite number of particles (without periodic boundaries).
<img src='http://www.etano.net/phys466/Grr-beads_021_m32_grr.png' />