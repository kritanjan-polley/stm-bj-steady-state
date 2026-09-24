# STM-BJ Steady-State Simulations

This repository contains simulation codes for steady-state transport in a scanning tunneling microscope break-junction (STM-BJ) model coupled to vibrational and cavity degrees of freedom using a quasiclassical mapping Hamiltonian approach.

This Github repo currently includes three methods:

- `semiclassicalMapping`:  Trajectory-based semiclassical mapping dynamics.
- `lindbladMasterEq`: Lindblad master-equation calculations using `QuantumToolbox.jl`.
- `HEOMcalculations`:  Hierarchical equations of motion (HEOM) calculations using `HierarchicalEOM.jl`.

## Running the Mapping Hamiltonian Codes

```bash
cd semiclassicalMapping
gfortran -cpp -O2 -ffast-math -funroll-loops \
    params.f90 functions.f90 propagator.f90 \
    initial.f90 main.f90 memory.c -framework Accelerate
```


## Running the Lindblad Master-Equation Codes

```bash
cd lindbladMasterEq
julia masterSTMBj3.jl
julia masterBareBJ.jl
```

## Running the HEOM Code

```bash
cd HEOMcalculations
julia HEOMjl.jl
```
