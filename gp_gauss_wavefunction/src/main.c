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

	f64 lambda = -1.0;
	i64 N  = 4;
	f64 g0 = lambda/((f64)N-1.0);

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

        .num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.5,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
	struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);

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
				  default_out*default_out,
				  random_out*random_out);
		}

		fclose(fd);
	}

	sbmf_shutdown();
}
