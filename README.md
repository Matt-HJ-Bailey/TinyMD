# TinyMD
A Molecular Dynamics simulation program to study Lennard-Jones fluids

## Features
Sets up an initial cubic lattiace in a simulation box
Implements a few different thermostatting algorithms:
* Berendsen
* Andersen
* Velocity Rescaling

Reads key:value pairs from a file titled `input.inpt`.
The following keys are read:
### Lennard-Jones Coefficients
* `LJ_sigma`
* `LJ_epsilon`
* `LJ_cutoff`
### Box dimensions
* `length_x`
* `length_y`
* `length_z`
* `num_x`
* `num_y`
* `num_z`
### Simulation parameters
* `mass`
* `temperature`
* `timesep`
* `num_rescale`
### Output parameters
* `num_output`
* `num_steps`
* `position_file`

## Requirements
* Eigen
* C++17
* CMake
