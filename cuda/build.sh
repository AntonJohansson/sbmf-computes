#!/bin/sh

CUDADIR=/opt/cuda/

gcc -o sum sum.c
${CUDADIR}/bin/nvcc -c -Xcompiler -c -o sumcu.o sum.cu -I ${CUDADIR}/include
