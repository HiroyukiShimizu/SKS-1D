# SKS-1D
SKS-1D is a Fortran90 code designed to solve the equations for a 1-D two-layer depth-averaged model of pyroclastic density currents [1, 2].


### Authors
Hiroyuki A. Shimizu (NIED), Takehiro Koyaguchi (ERI, UTokyo), Yujiro J. Suzuki (ERI, UTokyo)


### Installation and Test
(1) Installing gfortran, lapack, blas, and gnuplot

(2) Compiling SKS-1D by Makefile

    $ make

(3) Test

    $ cd Benchmark_dilute-PDC
    $ sh Allrun.sh


### Acknowledgements
The development of SKS-1D:1.0 has been partially funded by KAKENHI Grant Numbers JP17H02949 and JP21K14018 and by MEXT's "Integrated Program for Next Generation Volcano Research and Human Resource Development."


### References
[1] Shimizu H.A., Koyaguchi T., Suzuki Y.J. (2019) The run-out distance of large-scale pyroclastic density currents: A two-layer depth-averaged model. J. Volcanol. Geotherm. Res., 381, 168-184. https://doi.org/10.1016/j.jvolgeores.2019.03.013

[2] Shimizu H.A., Koyaguchi T., Suzuki Y.J., Brosch E., Lube G., Cerminara M. (2021) Validation of a two-layer depth-averaged model by comparison with an experimental dilute stratified pyroclastic density current. Bull. Volcanol., 83, 73. https://doi.org/10.1007/s00445-021-01493-w
