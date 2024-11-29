#!/bin/sh
#PJM -L rscgrp=regular-o
#PJM -L node=2x2x2
#PJM --mpi proc=1
#PJM --omp thread=16
#PJM -L elapse=72:00:00
#PJM -g gd22
#PJM -j

mpiexec ./a.out

