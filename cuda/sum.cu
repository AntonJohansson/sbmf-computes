#include <stdlib.h>
#include <assert.h>

#define N 100000000

__global__ void add(float* c, float* a, float* b, unsigned int n) {
	for (unsigned int i = 0; i < n; ++i) {
		c[i] = a[i] + b[i];
	}
}

int main() {
	float* a = (float*)malloc(N*sizeof(float));
	float* b = (float*)malloc(N*sizeof(float));
	float* c = (float*)malloc(N*sizeof(float));

	for (unsigned int i = 0; i < N; ++i) {
		a[i] = 1.0;
		b[i] = 2.0;
	}

	float* da;
	float* db;
	float* dc;
	cudaMalloc(&da, N*sizeof(float));
	cudaMalloc(&db, N*sizeof(float));
	cudaMalloc(&dc, N*sizeof(float));

	cudaMemcpy(da, a, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(db, b, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dc, c, N*sizeof(float), cudaMemcpyHostToDevice);

	add<<<1,1>>>(dc,db,da,N);

	cudaMemcpy(a, da, N*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(b, db, N*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(c, dc, N*sizeof(float), cudaMemcpyDeviceToHost);

	for (unsigned int i = 0; i < N; ++i) {
		assert(c[i] == 3.0f);
	}

	return 0;
}
