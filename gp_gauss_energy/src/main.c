#include <sbmf/sbmf.h>
#include <math.h>

#include <stdio.h>

#define PERTURBATION(x) 2*gaussian(x, 0, 0.1)
//#define PERTURBATION(x) (1.0/(cosh(2.0*x)*cosh(2.0*x)))

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

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_set_thread_storage_size(256*1024*1024);
	//OMEGA = 0.1;

	struct nlse_guess random_guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
	struct nlse_guess* default_guesses = NULL;

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 15000,
		.max_quadgk_iters = 500,
		//.abs_error_tol = 1e-9,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.75,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	//f64 lambda = -1.0;
	i64 Ns[] = {67,70,72,75};
	//i64 Ns[] = {127,132,137};
	f64 gs[] = {-0.01};

	for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {

		i64 N = Ns[i];
		for (u32 j = 0; j < sizeof(gs)/sizeof(gs[0]); ++j) {
			f64 g0 = gs[j];

			FILE* fd = fopen("out", "a");
			fprintf(fd, "%ld\t", N);

			sbmf_init();
			{
				struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
				if (!gp_default_res.converged)
					break;
				f64 Egp_default = grosspitaevskii_energy(settings, gp_default_res.coeff_count, component_count, gp_default_res.coeff, &N, &g0);
				f64 bmf_default = 0;//= bestmf_find_fractional_occupation(settings, N, g0, default_guesses);
				struct pt_result rs_default_ptres = rspt_1comp_cuda_new(&settings, gp_default_res, 0, g0, N);
				struct pt_result en_default_ptres = enpt_1comp_cuda_new(&settings, gp_default_res, 0, g0, N);
				fprintf(fd, "%.10f\t%.10lf\t", Egp_default, bmf_default);
				fprintf(fd, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t",
						rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2,
						rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2 + rs_default_ptres.E3,
						en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2,
						en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2 + en_default_ptres.E3
						);
			}
			sbmf_shutdown();

			sbmf_init();
			{
				struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);
				if (!gp_random_res.converged)
					break;
				f64 Egp_random  = grosspitaevskii_energy(settings, gp_random_res.coeff_count, component_count, gp_random_res.coeff, &N, &g0);
				f64 bmf_random = 0;//bestmf_find_fractional_occupation(settings, N, g0, random_guesses);
				struct pt_result rs_random_ptres = rspt_1comp_cuda_new(&settings, gp_random_res, 0, g0, N);
				struct pt_result en_random_ptres = enpt_1comp_cuda_new(&settings, gp_random_res, 0, g0, N);
				fprintf(fd, "%.10f\t%.10lf\t", Egp_random, bmf_random);
				fprintf(fd, "%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2,
						rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2 + rs_random_ptres.E3,
						en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2,
						en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2 + en_random_ptres.E3
						);
			}
			sbmf_shutdown();

			fclose(fd);
		}
	}
}
