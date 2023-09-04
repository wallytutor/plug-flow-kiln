# Rotary Kiln Model

Implementation of rotary kiln model for simulation of thermal processing of materials.

Details of available resources and a user guide are provided in the [documentation](docs/).

## Model validation

Partial model validation is provided under [`validation/`](validation) directory. For now it covers gas properties and mass action kinetics computed with help of Cantera and the proper functioning of materials constructors and radiation model. Further work is required for full coverage and physical correcteness check.

## Governance and development

Currently the code is written in Octave using package [CasADi](https://web.casadi.org/) for handling symbolics and solution of nonlinear equations. The goal of Octave was initially to allow fast prototyping in an environment with built-in support of vectorized operations and numerical computations. Although Octave is pretty mature and convenient for mathematical modeling, its I/O and graphical features are poor for production-grade code. Also its classes have limited features. It was thus decided to adopt the following approach:

1. All new features are added to the Octave code base until a Python translation is fully implemented.

1. Python implementation will also rely on [CasADi](https://web.casadi.org/) and just a few major well supported packages.

1. Once a Python code base is established, Octave support will be dropped and code placed in a legacy directory.

1. Python code that is considered stable will be progressively migrated to C++ and a Python API will be added.

1. In the long-term Python will be the official scripting language for application programs with code base fully in C++.

1. All code must be Linux (gcc) and Windows (MingW64) compatible.

1. Long-term goal is to implement both a 2-D OpenFOAM solver based on this work and a set of Ansys Fluent UDF's for rotary kiln simulation.
