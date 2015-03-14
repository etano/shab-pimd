# Thermostats #

As written in the introduction, our molecular dynamics simulation is in the constant energy micro-canonical ensemble. To sample the partition function at a constant temperature canonical ensemble a thermostat is usually required. Note that necessarily this is unphysical since one is adding degrees of freedom to the system which are not truly present and changing the conserved quantity (the Hamiltonian). Nevertheless, there has been much success with correctly reproducing the correct canonical distribution through thermostatting. We compare two methods:
  * The [Nose\_Hoover\_Thermostat](Nose_Hoover_Thermostat.md) deterministically introduces new degrees of freedom that lead to a constant temperature. It is considered the _gold standard_ of thermostatting.
  * The [Langevin\_Thermostat](Langevin_Thermostat.md) makes use of the Langevin Equation to mimic a heat bath and it thus stochastic. It employs the Einstein relation between friction and the random force to set the temperature.

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)