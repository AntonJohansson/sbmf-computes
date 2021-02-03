#!/bin/sh

CUDADIR=/opt/cuda/

gcc -o sum sum.c
${CUDADIR}/bin/nvcc -o sumcu sum.cu -I ${CUDADIR}/include
