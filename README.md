# STM-BJ Steady-State Simulations

This repository contains simulation codes for steady-state transport in a
scanning tunneling microscope break-junction (STM-BJ) model coupled to
vibrational and cavity degrees of freedom using a quasiclassical mapping Hamiltonian 
approach.

The codebase currently includes three complementary approaches:

- `semiclassicalMapping/`: Fortran trajectory-based semiclassical mapping
  dynamics.
- `lindbladMasterEq/`: Julia Lindblad master-equation calculations using
  `QuantumToolbox.jl`.
- `HEOMcalculations/`: Julia hierarchical equations of motion (HEOM)
  calculations using `HierarchicalEOM.jl`.

## Running the Mapping Hamiltonian Codes

```bash
cd semiclassicalMapping
gfortran -cpp -O2 -ffast-math -funroll-loops \
    params.f90 functions.f90 propagator.f90 \
    initial.f90 main.f90 ../lindbladMasterEq/memory.c -framework Accelarate
```


## Running the Lindblad Master-Equation Codes

```bash
cd lindbladMasterEq
cc -dynamiclib -fPIC memory.c -o libmemory.dylib  # macOS
cc -shared -fPIC memory.c -o libmemory.so # linux
```
