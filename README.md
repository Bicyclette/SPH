# SPH project

The idea of this project is to compare the speed of different SPH implementations.

## HOWTO

Be sure to have the following dependencies installed in your machine:
- OpenMP (for parallel execution)

To run the project, simply run the following lines in a terminal:

- Linux
```sh
cmake -B build -S .
cmake --build build --config Release
cd build
./sph
```
