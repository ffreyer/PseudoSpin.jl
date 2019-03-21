# PseudoSpin

This is the Monte-Carlo simulation behind the paper
[Two-stage multipolar ordering in PrT<sub>2</sub>Al<sub>20</sub> Kondo materials](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.115111). It features an effective classical model on the diamond lattice with

* J<sub>1</sub>: nearest neighbor coupling (xxz anisotropy)
* J<sub>2</sub>: second nearets neighbor coupling (xxz anisotropy)
* J<sub>3</sub>: third nearest neighbor coupling (xxz anisotropy)
* Γ: three-spin chain interaction
* K: four-spin chain interaction
* ζ: spatially anisotropic second neighbor coupling
* h: magnetic field

At the time of the first paper J<sub>1</sub>, J<sub>2</sub>, K and h were used.

---

### Basic usuage

The front-end function to run the simulation is `simulate!()`. It takes a large number of keyword arguments to set up the details of a simulation. The following **keyword arguments** are available

###### Output file directory:

The full path will `path * folder * filename * ".part"`. `folder` is often useful to seperate results by parameters when running a large number of simulations.
* `path = ""`: path to output file
* `folder = ""`: name of output folder
* `filename = "T"`: filename

###### Lattice setup:

The lattice generation currently performs a spatial search to find neighbors, thus the first parameter. The path generation is scattered, thus only the four-spin interaction is mentioned here.

* `neighbor_search_depth = 3`: Number of neighbor orders used (e.g. three for third nearest neighbors)
* `do_paths = true`: whether to include 4-spin paths
* `L = 6`: linear system size (resulting in 2L<sup>3</sup> sites)
* `spins = rand_spin(2*L^3)`: initial vector of spins

###### Simulation parameters:

Note that some keyword arguments take priority over other. For example, `J1` and `lambda` combine to `J1s`, which is then used in the simulation. If `J1s` is given, `J1` (and `lambda`) is ignored.

* `T = 1.0`: Temperature for serial runs
* `Ts = [T]`: Temperatures for parallel tempering runs (priority over T)
* `J1 = 0.`: (isotropic) nearest neighbor coupling
* `J2 = 0.`: (isotropic) second nearest neighbor coupling
* `J3 = 0.`: (isotropic) third nearest neighbor coupling
* `lambda = 0.`: joint anisotropy for J1, J2, J3
* `J1s = (J1, lambda*J1)`: (anisotropic) nearest neighbor coupling (xy, z) (priority over J1)
* `J2s = (J2, lambda*J2)`: (anisotropic) second nearest neighbor coupling (xy, z) (priority over J2)
* `J3s = (J3, lambda*J3)`: (anisotropic) third nearest neighbor coupling (xy, z) (priority over J3)
* `K = 0.`: four-spin chain coupling
* `h = Point3(0.)`: magnetic field strength
* `g = 0.`: three-spin chain coupling
* `zeta = 0.`: spatially anisotropic second nearest neighbor coupling

###### Thermalization - Temperature generator

The `Tgen_method` (temperature generation method) is an iterator returning temperatures used in the thermalization step of the simulation. The options here include `ConstantT`, which continously returns `T` given above, and `Freezer` which performs simulated annealing. More specifically `Freezer` generated temperatures following an exponential decay, modulated with a sine, starting at `Freeze_temperature` and becoming constant after `N_switch` sweeps at `T`. (Both work with serial and parallel simulations)

* `TH_sweeps = 2_000_000`: Number of thermalization sweeps
* `TGen_method = Freezer`: Temperature generator for thermalization. (ConstantT or Freezer, see above)
* `N_switch = div(TH_sweeps, 2)`: Time scale for exponential decay of `Freezer` (ignored for `ConstantT`)
* `Freeze_temperature = 1.5*maximum(Ts)`: Initial temperature for `Freezer` (ignored for `ConstantT`)

###### Thermalization - Parallel Tempering

* `batch_size = 10`: Number of sweeps between parallel tempering steps
* `thermalizer_method = NoParallelTempering`: Thermalization method used (`NoParallelTempering` or `ParallelTempering`)


###### Measurement
* `ME_sweeps = 5_000_000`: Number of sweeps for the measurement
