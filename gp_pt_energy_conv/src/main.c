#include <sbmf/sbmf.h>

#include <stdio.h>

#define NA 4
#define NB 0

#define GAA (1.0/3.0)
//#define GAA (-2.0/((f64)NA-1))
#define GAB (+1.0/((f64)NB))
#define GBA (+1.0/((f64)NA))
#define GBB (-2.0/((f64)NB-1))

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
#else
	struct nlse_guess* guesses = NULL;
#endif

	f64 g0[] = {
		GAA, GAB,
		GBA, GBB
	};

	i64 occupations[] = {NA,NB};

	f64 g0s[] = {-1.0/6.0, -1.0/3.0, 1.0/3.0, 1.0/6.0};
	u32 bs[] = {4,8,12,16,24,32,48,64,80};
	struct nlse_settings settings = {
        //.spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.error_tol = 1e-14,

		.num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(g0s)/sizeof(g0s[0]); ++j) {
		g0[0] = g0s[j];

		for (u32 i = 0; i < sizeof(bs)/sizeof(bs[0]); ++i) {
			u32 b = bs[i];
			settings.num_basis_funcs = b;

			sbmf_init();
			struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
			f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);

			struct pt_result rs_ptres = rspt_1comp_cuda_new(&settings, res, 0, g0[0], occupations[0]);
			struct pt_result en_ptres = enpt_1comp_cuda_new(&settings, res, 0, g0[0], occupations[0]);

			//printf("E0:          %.15lf\n", ptres.E0);
			//printf("E1:          %.15lf\n", ptres.E1);
			//printf("E2:          %.15lf\n", ptres.E2);
			//printf("E3:          %.15lf\n", ptres.E3);
			//printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
			//printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
			//printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

			{
				char buf[256];
				snprintf(buf, 256, "out_g%.2lf", g0[0]);

				FILE* fd = fopen(buf, "a");
				fprintf(fd, "%u\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						b,
						Egp,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3,
						en_ptres.E0+en_ptres.E1+en_ptres.E2,
						en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3
					   );
				fclose(fd);
			}

			sbmf_shutdown();
		}
	}

}
