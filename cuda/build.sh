#!/bin/sh

gcc -o sum sum.c
/opt/cuda/bin/nvcc -o sumcu sum.cu
