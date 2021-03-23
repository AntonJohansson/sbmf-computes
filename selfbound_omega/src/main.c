#include <sbmf/sbmf.h>

#include <stdio.h>

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define USE_RANDOM_GUESS 0












void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i) {
		//out[i] = (1 - 0.5 * in[i]*in[i]);
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);// + gaussian(in[i] - 1.0, 0.0, 0.2);
	}
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}



int main() {
	sbmf_set_log_callback(log_callback);

#if USE_GAUSSIAN_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian0,
		},
		[1] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian1,
		},
	};
#elif USE_RANDOM_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
#else
	struct nlse_guess* guesses = NULL;
#endif

	struct nlse_settings settings = {
		.max_iterations = 2000,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 48,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.hamiltonian_mixing = 0.7,

		.orbital_choice = NLSE_ORBITAL_LOWEST_ENERGY,

		.gk=gk20
    };

	const u32 component_count = 2;

	i64 N = 4;
	//f64 Os[] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
	f64 Os[] = {0.005};
	f64 fs[] = {-4,-5,-6,-7,-8,-9,-10};

	f64 lambda = 0.5;
	for (u32 j = 0; j < sizeof(fs)/sizeof(fs[0]); ++j) {
		f64 factor = fs[j];

		for (u32 i = 0; i < sizeof(Os)/sizeof(Os[0]); ++i) {
			OMEGA = Os[i];
			i64 occupations[] = {N,N};
			f64 g0[] = {
				lambda/((f64)N-1), factor*lambda/((f64)N-1),
				factor*lambda/((f64)N-1),  lambda/((f64)N-1)
			};

			sbmf_init();
			struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
			if (!res.converged)
				break;

			f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
			printf("\nfull energy: %lf\n", Efull);

			// out Energy
			{
				char buf[256];
				snprintf(buf, 256, "out_f%lf_E", fs[j]);
				FILE* fd = fopen(buf, "a");
				fprintf(fd, "%lf\t%.10e\n",
						OMEGA,
						Efull
						);
				fclose(fd);
			}

			// out wf
			{
				char buf[256];
				snprintf(buf, 256, "out_f%lf_W_w%lf", fs[j], Os[i]);
				FILE* fd = fopen(buf, "a");
				for (f64 x = -10; x < 10; x += 0.01) {
					f64 outA, outB;
					ho_sample(res.coeff_count, &res.coeff[0], 							1, &outA, &x);
					ho_sample(res.coeff_count, &res.coeff[res.coeff_count], 1, &outB, &x);
					fprintf(fd, "%lf\t%lf\t%lf\n", x, outA*outA, outB*outB);
				}

				fclose(fd);
			}

			//struct pt_result rspt = rspt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
			//struct pt_result enpt = enpt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
			//{
			//	FILE* fd = fopen("out", "a");
			//	fprintf(fd, "%lf\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
			//			OMEGA,
			//			Efull,
			//			rspt.E0+rspt.E1+rspt.E2,
			//			rspt.E0+rspt.E1+rspt.E2+rspt.E3,
			//			enpt.E0+enpt.E1+enpt.E2,
			//			enpt.E0+enpt.E1+enpt.E2+enpt.E3
			//			);
			//	fclose(fd);
			//}

			sbmf_shutdown();
		}
	}



}
