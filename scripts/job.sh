#!/bin/sh
#PJM -L rscgrp=regular-o
#PJM -L node=1x1x1
#PJM --mpi proc=1
#PJM --omp thread=16
#PJM -L elapse=12:00:00
#PJM -g gn38
#PJM -j

mpiexec ./a.out

