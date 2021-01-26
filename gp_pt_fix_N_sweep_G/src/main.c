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
	struct nlse_guess* default_guesses = NULL;

	i64 N = 100;
	f64 l0s[] = {-5,-4,-3,-2,0,-1,1,2,3,4,5};
	//f64 l0s[] = {};
	//i64 range_N[] = {100, 500, 15000};
	struct nlse_settings settings = {
        //.spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_integration_evals = 1e5,
		.error_tol = 1e-14,

        .num_basis_funcs = 16,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.5,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(l0s)/sizeof(l0s[0]); ++j) {
		sbmf_init();

		struct nlse_result res;
		struct bestmf_result bmf_gaussian_res;
		struct bestmf_result bmf_default_res;

		f64 g0 = l0s[j]/((f64)N-1.0);

		res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
		f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

		bmf_gaussian_res = best_meanfield(settings, N, g0, gaussian_guesses);
		bmf_default_res = best_meanfield(settings, N, g0, default_guesses);

		struct pt_result rs_ptres = rayleigh_schroedinger_pt_rf(settings, res, 0, &g0, &N);
		struct pt_result en_ptres = en_pt_rf(settings, res, 0, &g0, &N);

		//printf("E0:          %.15lf\n", ptres.E0);
		//printf("E1:          %.15lf\n", ptres.E1);
		//printf("E2:          %.15lf\n", ptres.E2);
		//printf("E3:          %.15lf\n", ptres.E3);
		//printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		//printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		//printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		{
			char buf[256];
			snprintf(buf, 256, "out", l0s[j]);

			FILE* fd = fopen(buf, "a");
			fprintf(fd, "%.10lf\t%.10lf\t%.10f\t%.10f\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
					l0s[j],
					Egp,
					bmf_gaussian_res.energy								- Egp,
					bmf_default_res.energy 								- Egp,
					rs_ptres.E0+rs_ptres.E1+rs_ptres.E2					- Egp,
					rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3		- Egp,
					en_ptres.E0+en_ptres.E1+en_ptres.E2					- Egp,
					en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3		- Egp
				   );
			fclose(fd);
		}


		sbmf_shutdown();
	}

}
