## HanleRT-TIC

HanleRT-TIC is a numerical code for solving the problem of the generation and
transfer of polarized radiation in one-dimensional (1D and 1.5D) model
atmospheres for arbitrary atomic models. HanleRT-TIC also allows for the
solution of the inversion problem with its Tenerife Inversion Code (TIC)
module.

The code can take into account atomic and scattering polarization, quantum
interference among energy levels pertaining to a given atomic term, partial
frequency redistribution effects due to the photon coherence in scattering
processes, the impact of dynamics, and the impact of magnetic fields via the
joint action of the Hanle and Zeeman effects (incomplete Paschen-Back regime).

Apart from solving the radiation transfer problem in optically thick
plane-parallel model atmospheres, HanleRT-TIC allows for the solution of the
coronal line emission (CLE) problem, in which the radiation field is assumed
to be dominated by the underlying stellar disk.

HanleRT-TIC is written in standard Fortran 2008, parallelized with the
[OpenMPI](https://www.open-mpi.org/) and [OpenMP](https://www.openmp.org/)
libraries. The code also has some parsing routines written in python to
allow for more flexible input formats.

You can find the documentation [here](https://tdpa.gitlab.io/hanlert-tic/)

## Terms of use

You can use HanleRT-TIC freely in your scientific research.

In case you publicly present your results obtained with the help of the
of the synthesis capability of the code (HanleRT), we ask you to quote the paper
[del Pino Alemán et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...830L..24D/abstract).
If using the inversion capabilities (HanleRT-TIC), we ask you to also quote
[Li et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...933..145L/abstract).
In both cases, please, include a link to the code repository
[https://gitlab.com/TdPA/hanlert-tic](https://gitlab.com/TdPA/hanlert-tic).

Authors
=======

* Tanausú del Pino Alemán

* Hao Li

* Roberto Casini

## License

MIT License

Copyright (c) 2023 Tanausú del Pino Alemán

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

