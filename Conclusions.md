# Conclusions #
We are able to use PIMD to model a quantum mechanical system. At high temperatures, PIMD simulations approach a classical result where all of the beads move close together. At lower temperatures, the beads spread out, corresponding to a quantum mechanical regime. Our simulations clearly showed the relationship between bead spreading and classical or quantum mechanical behavior.

We found that the number of time slices, or beads on the necklace, required for quantum convergence varies depending on the system. Also, the time step for high temperatures needs to be smaller than what we used to run our simulations. However, the focus of our project was on the comparison of different method combinations, so we used the same time step in every simulation.

The best combination of methods is clearly the Langevin thermostat with the staging transformation. This conclusion is based upon the time that it took to run the simulation and the errors in the values predicted by the simulation. The next best combination is the Langevin thermostat with the normal mode transformation. Even with our optimizations, the Nose-Hoover thermostat in combination with any of the transformations did not compare to the use of the Langevin thermostat with the staging transformation.

# Future Directions #

In the future, it would be interesting to apply the PIMD approach to an actual system using our optimal settings. Analysis of systems such as the interaction of hydrogen with a metal lattice would benefit from the use of the PIMD method, because it would incorporate quantum mechanical effects due to the small mass of the hydrogen atom. For this application, PIMD would be combined with an ab initio method such as DFT. For futher comparisons, PIMD should be compared with other methods such as the path integral monte carlo method, looking at efficiency.


---

**NEXT SECTION**: [References](References.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)