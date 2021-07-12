# SKS-1D:1.0

SKS-1D:1.0 is a Fortran90 code designed to solve the equations for a 1-D two-layer depth-averaged model of pyroclastic density currents [1][2].

### Authors
Hiroyuki A. Shimizu (NIED), Takehiro Koyaguchi (ERI, UTokyo), Yujiro J. Suzuki (ERI, UTokyo)

### Installation and Test
(1) Installing gfortran, lapack, blas, and gnuplot
(2) Compiling SKS-1D:1.0 by Makefile
    $ make
(3) Test
    $ cd Benchmark
    $ sh Allrun.sh

### Acknowledgements
The development of SKS-1D:1.0 has been partially funded by KAKENHI Grant Numbers JP17H02949 and JP21K14018 and by MEXT's "Integrated Program for Next Generation Volcano Research and Human Resource Development."

### References
[1] Shimizu H.A., Koyaguchi T., Suzuki Y.J. (2019) The run-out distance of large-scale pyroclastic density currents: A two-layer depth-averaged model. J. Volcanol. Geotherm. Res., 381, 168-184.
[2] Shimizu H.A., Koyaguchi T., Suzuki Y.J., Brosch E., Lube G., Cerminara M. (in submit) Validation of a two-layer depth-averaged model by comparison with an experimental dilute stratified pyroclastic density current. Bull. Volcanol.
