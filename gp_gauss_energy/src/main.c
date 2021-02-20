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

struct full_energy_integrand_params {
    u32 coeff_count;
    f64* coeff_a;
    f64* coeff_b;
    struct basis basis;
};

/* Integrand of the form |a|^2|b|^2 */
void norm_overlap_integrand(f64* out, f64* in, u32 len, void* data) {
    struct full_energy_integrand_params* p = data;

    f64 sample_a[len];
    p->basis.sample(p->coeff_count, p->coeff_a, len, sample_a, in);

    f64 sample_b[len];
    p->basis.sample(p->coeff_count, p->coeff_b, len, sample_b, in);

    for (u32 i = 0; i < len; ++i) {
        out[i] = sample_a[i]*sample_a[i]*sample_b[i]*sample_b[i];
    }
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
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 32,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.5,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	f64 lambda = -0.0;
	//i64 Ns[] = {4,6,8,12,16,20,24,28,32};
	i64 Ns[] = {8};

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

		{
			struct full_energy_integrand_params p = {
				.coeff_count = gp_default_res.coeff_count,
				.coeff_a = &gp_default_res.coeff[0],
				.coeff_b = &gp_default_res.coeff[gp_default_res.coeff_count],
				.basis = settings.basis,
			};

			struct quadgk_settings int_settings = {
				.max_iters = 500,
				.abs_error_tol = 1e-15,
				.userdata = &p,
				.gk = gk20
			};

			u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];

			struct quadgk_result ires;
			quadgk_infinite_interval(norm_overlap_integrand, &int_settings, quadgk_memory, &ires);
			f64 B = ires.integral;

			p.coeff_a = &gp_default_res.coeff[0];
			quadgk_infinite_interval(norm_overlap_integrand, &int_settings, quadgk_memory, &ires);
			f64 A0 = ires.integral;

			//f64 NC = (Egp_default - Egp_random)/(3.0*B - A0);
			f64 NC = (gp_default_res.energy[1] - gp_default_res.energy[0])/(3.0*B - A0);
			printf("..... NC: %lf\n", NC);
		}


#if 0
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
#endif

		sbmf_shutdown();
	}
	fclose(default_fd);
	fclose(random_fd);
}
