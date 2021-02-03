#include <stdlib.h>
#include <assert.h>

#define N 100000000

void main() {
	float* a = malloc(N*sizeof(float));
	float* b = malloc(N*sizeof(float));
	float* c = malloc(N*sizeof(float));

	for (unsigned int i = 0; i < N; ++i) {
		a[i] = 1.0;
		b[i] = 2.0;
	}

	for (unsigned int i = 0; i < N; ++i) {
		c[i] = a[i] + b[i];
	}

	for (unsigned int i = 0; i < N; ++i) {
		assert(c[i] == 3.0f);
	}
}
