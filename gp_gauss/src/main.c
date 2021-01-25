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
	sbmf_init();

	struct nlse_guess random_guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
	struct nlse_guess* default_guesses = NULL;

	f64 g0 = 1.0/3.0;
	i64 N  = 4;

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_integration_evals = 1e5,
		.error_tol = 1e-14,

        .num_basis_funcs = 64,
		.basis = ho_basis,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
	struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);

	f64 Egp_default = grosspitaevskii_energy(settings, gp_default_res.coeff_count, component_count, gp_default_res.coeff, &N, &g0);
	f64 Egp_random  = grosspitaevskii_energy(settings, gp_random_res.coeff_count, component_count, gp_random_res.coeff, &N, &g0);

	struct bestmf_result bmf_default = best_meanfield(settings, N, g0, default_guesses);
	struct bestmf_result bmf_random = best_meanfield(settings, N, g0, random_guesses);

	struct pt_result rs_default_ptres = rayleigh_schroedinger_pt_rf(settings, gp_default_res, 0, &g0, &N);
	struct pt_result en_default_ptres = en_pt_rf(settings, gp_default_res, 0, &g0, &N);

	struct pt_result rs_random_ptres = rayleigh_schroedinger_pt_rf(settings, gp_random_res, 0, &g0, &N);
	struct pt_result en_random_ptres = en_pt_rf(settings, gp_random_res, 0, &g0, &N);

	{
		FILE* fd = fopen("out_E", "a");
		fprintf(fd, "# mu\tE\tBMF\tRS2\tRS3\tEN2\tEN3\n");
		fprintf(fd, "#default\n%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
				gp_default_res.energy[0],
				Egp_default,
				bmf_default.energy,
				rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2,
				rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2 + rs_default_ptres.E3,
				en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2,
				en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2 + en_default_ptres.E3
				);
		fprintf(fd, "#random\n%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
				gp_random_res.energy[0],
				Egp_random,
				bmf_random.energy,
				rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2,
				rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2 + rs_random_ptres.E3,
				en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2,
				en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2 + en_random_ptres.E3
				);

		fclose(fd);
	}

	{
		FILE* fd = fopen("out_W", "a");
		fprintf(fd, "# x\tpotential\tdefault\trandom\n");

		for (f64 x = -3; x <= 3; x += 0.01) {
		  /* Sampling one x at a time is inefficient, but it doesnt really matter here */
		  f64 default_out = 0, random_out = 0;
		  ho_sample(gp_default_res.coeff_count, &gp_default_res.coeff[0], 1, &default_out, &x);
		  ho_sample(gp_random_res.coeff_count, &gp_random_res.coeff[0], 1, &random_out, &x);
		  f64 pot = 0.5*x*x + PERTURBATION(x);
		  fprintf(fd, "%lf\t%lf\t%lf\t%lf\n", x, pot,
				  gp_default_res.energy[0] + default_out*default_out,
				  gp_random_res.energy[0] + random_out*random_out);
		}

		fclose(fd);
	}

	sbmf_shutdown();
}
