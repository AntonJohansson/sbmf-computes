#include <sbmf/sbmf.h>

#include <stdio.h>

#define USE_GAUSSIAN_GUESS 0

//#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)
#define PERTURBATION(x) 0.0
//#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = PERTURBATION(in_x[i]);
    }
}





















void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}





void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}




int main() {
	sbmf_set_log_callback(log_callback);

	struct nlse_guess gaussian_guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian0,
		},
		[1] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian1,
		},
	};
	struct nlse_guess random_guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
	struct nlse_guess* default_guesses = NULL;

	f64 l0s[] = {-0.5, -1.0, 1.0, 0.5};
	//f64 l0s[] = {};
	//i64 range_N[] = {100, 500, 15000};
	//i64 Ns[] = {4,6};
	i64 Ns[] = {
		4,
		6,
		10,
		20,
		30,
		40,
		60,
		80,
		100,
		120,
		140,
		160
	};
	//i64 Ns[] = {4,30,100};
	struct nlse_settings settings = {
		//.spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_quadgk_iters = 1000,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(l0s)/sizeof(l0s[0]); ++j) {

		for (u32 k = 0; k < sizeof(Ns)/sizeof(Ns[0]); ++k) {
			sbmf_init();

			i64 N = Ns[k];
			f64 g0 = l0s[j]/((f64)N-1.0);


			struct nlse_result res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
			f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

			//{
			//	f64 bmf_gaussian  = bestmf_find_fractional_occupation(settings, N, g0, gaussian_guesses);
 			//	f64 bmf_default   = bestmf_find_fractional_occupation(settings, N, g0, default_guesses);
 			//	f64 bmf_random    = bestmf_find_fractional_occupation(settings, N, g0, random_guesses);
			//	{
			//		char buf[256];
			//		snprintf(buf, 256, "out_bmf_l%.2lf", l0s[j]);

			//		FILE* fd = fopen(buf, "a");
			//		fprintf(fd, "%ld\t%.10e\t%.10e\t%.10e\n",
			//				N,
			//				bmf_default,
			//				bmf_gaussian,
			//				bmf_random
			//			   );
			//		fclose(fd);
			//	}
			//}

			struct pt_result rs_ptres = rspt_1comp_cuda_new(&settings, res, 0, g0, N);
			struct pt_result en_ptres = enpt_1comp_cuda_new(&settings, res, 0, g0, N);

			{
				char buf[256];
				snprintf(buf, 256, "out_l%.2lf", l0s[j]);

				FILE* fd = fopen(buf, "a");
				fprintf(fd, "%ld\t%.10f\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						N,
						Egp,
						rs_ptres.E0+rs_ptres.E1,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3,
						en_ptres.E0+en_ptres.E1,
						en_ptres.E0+en_ptres.E1+en_ptres.E2,
						en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3
					   );
				fclose(fd);
			}

			sbmf_shutdown();
		}

	}

}
