# nbody_jl
This is a simple gravitational n-body simulation routine written in Julia. It uses Velocity Verlet integration to compute trajectories of n bodies with varying masses m. It outputs the body positions to a csv file (can be opened in almost all text editors or spreadsheet applications) and saves an svg image file of plotted trajectories. Simulations are set up with a simple text configuration file called (an example "input.txt" is provided).

Example plot outputs:

2-body example
![2_body_example](https://raw.githubusercontent.com/ilovematter/nbody_jl/f9fa3e6dbd704036233743ce21dbdbd4ecbace09/2_body_example.svg)

3-body example
![3_body_example](https://raw.githubusercontent.com/ilovematter/nbody_jl/5acf4d2681fd43893df5a9332697c28595f296b7/3_body_example.svg)
