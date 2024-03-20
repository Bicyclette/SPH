# SPH project

The idea of this project is to compare the speed of different SPH implementations.

## HOWTO

Be sure to have the following dependencies installed in your machine:
- OpenMP (for parallel execution)

To run the project, simply run the following lines in a terminal:

- Linux (for Arch users, first uncomment the lines in the `CMakeLists.txt`)
```sh
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release
./build/sph
```