#include <sbmf/sbmf.h>

#include <stdio.h>

#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)

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
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.5,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	f64 lambda = -1.0;
	i64 Ns[] = {4,6,8,12,16,20,24,28,32};

	FILE* default_fd = fopen("out_default_E", "a");
	FILE* random_fd = fopen("out_random_E", "a");

	fprintf(default_fd, "#N\tmu\tE\tBMF\tRS2\tRS3\tEN2\tEN3\n");
	fprintf(random_fd, "#N\tmu\tE\tBMF\tRS2\tRS3\tEN2\tEN3\n");

	for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {
		sbmf_init();

		i64 N = Ns[i];
		f64 g0 = lambda/((f64)N-1.0);

		struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
		struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);

		f64 Egp_default = grosspitaevskii_energy(settings, gp_default_res.coeff_count, component_count, gp_default_res.coeff, &N, &g0);
		f64 Egp_random  = grosspitaevskii_energy(settings, gp_random_res.coeff_count, component_count, gp_random_res.coeff, &N, &g0);

		struct bestmf_result bmf_default = best_meanfield(settings, N, g0, default_guesses);
		struct bestmf_result bmf_random = best_meanfield(settings, N, g0, random_guesses);

		struct pt_result rs_default_ptres = rspt_1comp_cuda_new(&settings, gp_default_res, 0, g0, N);
		struct pt_result en_default_ptres = enpt_1comp_cuda_new(&settings, gp_default_res, 0, g0, N);

		struct pt_result rs_random_ptres = rspt_1comp_cuda_new(&settings, gp_random_res, 0, g0, N);
		struct pt_result en_random_ptres = enpt_1comp_cuda_new(&settings, gp_random_res, 0, g0, N);

		{
			fprintf(default_fd, "%ld\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
					N,
					gp_default_res.energy[0],
					Egp_default,
					bmf_default.energy,
					rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2,
					rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2 + rs_default_ptres.E3,
					en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2,
					en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2 + en_default_ptres.E3
				   );
		}
		{
			fprintf(random_fd, "%ld\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
					N,
					gp_random_res.energy[0],
					Egp_random,
					bmf_random.energy,
					rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2,
					rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2 + rs_random_ptres.E3,
					en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2,
					en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2 + en_random_ptres.E3
				   );
		}

		sbmf_shutdown();
	}
	fclose(default_fd);
	fclose(random_fd);
}
