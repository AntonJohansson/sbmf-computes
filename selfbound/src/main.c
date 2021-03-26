#include <sbmf/sbmf.h>

#include <stdio.h>

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define USE_RANDOM_GUESS 1












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
		.max_iterations = 10000,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.hamiltonian_mixing = 0.9,

		.orbital_choice = NLSE_ORBITAL_LOWEST_ENERGY,

		.gk=gk20
    };

	OMEGA = 0.1;

	const u32 component_count = 2;

	i64 Ns[] = {4,6,10,15,20,40,60,80,100,120,140,160,180,200,300,400,500};
	//i64 Ns[] = {110};//,120,125,130,135,140,145,150,155,160,170,180,200,220,240,260,280,300};
	//i64 Ns[] = {200,225, 250, 275, 300, 325, 350, 375, 400,425, 450,475, 500};
	//i64 Ns[] = {200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000,3500,4000};
	//i64 Ns[] = {1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,4000,8000,16000};
	//i64 Ns[] = {2000,3000,4000,5000,6000,7000,8000,9000,10000};
	//f64 Os[] = {0.3, 0.2, 0.15, 0.1};
	//f64 Os[] = {0.1, 0.05, 0.025};
	f64 Os[] = {0.005};
	f64 gAB_factors[] = {-1.5};

	f64 lambda = 1.0;
	for (u32 k = 0; k < sizeof(gAB_factors)/sizeof(gAB_factors[0]); ++k) {
		f64 gAB_factor = gAB_factors[k];

		for (u32 j = 0; j < sizeof(Os)/sizeof(Os[0]); ++j) {
			OMEGA = Os[j];

			char buf[128];
			snprintf(buf, 128, "out_%lf_gab_%lf", Os[j], gAB_factor);
			{
				FILE* fd = fopen(buf, "a");
				fprintf(fd, "# N\tE\tRS2\tRS3\tEN2\tEN3\n");
				fclose(fd);
			}


			for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {
				i64 N = Ns[i];
				i64 occupations[] = {N,N};
				//f64 g0[] = {
				//	 lambda/((f64)N-1), gAB_factor*lambda/((f64)N-1),
				//	 gAB_factor*lambda/((f64)N-1),  lambda/((f64)N-1)
				//};

				const i64 max_N = 500;
				const f64 g = lambda/((f64)max_N - 1.0);
				f64 g0[] = {
					g, 						gAB_factor*g,
					gAB_factor*g, g
				};

				sbmf_init();
				struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
				if (!res.converged)
					break;
				f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
				printf("\nfull energy: %lf\n", Efull);

				struct pt_result rspt = rspt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
				struct pt_result enpt = enpt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
				{
					FILE* fd = fopen(buf, "a");
					fprintf(fd, "%ld\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
							N,
							Efull,
							rspt.E0+rspt.E1+rspt.E2,
							rspt.E0+rspt.E1+rspt.E2+rspt.E3,
							enpt.E0+enpt.E1+enpt.E2,
							enpt.E0+enpt.E1+enpt.E2+enpt.E3
							);
					fclose(fd);
				}

				sbmf_shutdown();
			}
		}
	}



}
