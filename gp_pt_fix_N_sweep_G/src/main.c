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
	//f64 l0s[] = {-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2};
	f64 l0s[] = {-4.5, -5};
	//f64 l0s[] = {};
	//i64 range_N[] = {100, 500, 15000};
	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.7,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(l0s)/sizeof(l0s[0]); ++j) {
		sbmf_init();

		struct nlse_result res;
		f64 bmf_gaussian_res;
		f64 bmf_default_res;

		f64 g0 = l0s[j]/((f64)N-1.0);

		res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
		f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

		bmf_gaussian_res = bestmf_find_fractional_occupation(settings, N, g0, gaussian_guesses);
		bmf_default_res = bestmf_find_fractional_occupation(settings, N, g0, default_guesses);

		struct pt_result rs_ptres = rspt_1comp_cuda_new(&settings, res, 0, g0, N);
		struct pt_result en_ptres = enpt_1comp_cuda_new(&settings, res, 0, g0, N);

		//printf("E0:          %.15lf\n", ptres.E0);
		//printf("E1:          %.15lf\n", ptres.E1);
		//printf("E2:          %.15lf\n", ptres.E2);
		//printf("E3:          %.15lf\n", ptres.E3);
		//printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		//printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		//printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		{
			FILE* fd = fopen("out", "a");
			fprintf(fd, "%.10lf\t%.10lf\t%.10f\t%.10f\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
					l0s[j],
					Egp,
					bmf_gaussian_res,
					bmf_default_res,
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
