#include <sbmf/sbmf.h>
#include <stdio.h>

#define NA 4

#define ORDER (10.0)

#define GAA (1.0/ORDER)

#define USE_GAUSSIAN_GUESS 0
#define USE_TF_GUESS 0

static f64 g0[] = { GAA };
static i64 occupations[] = {NA};

void perturbation_zero(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = 0.0;
    }
}

void perturbation_gauss(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = 2*gaussian(in_x[i], 0, 0.2);
    }
}

void tf(f64* out, f64* in, u32 len, void* data) {
	f64 mu = 0.5 * sqrt(3*g0[0]*(occupations[0]-1)/2);
	for (u32 i = 0; i < len; ++i) {
		out[i] = (mu - 0.5*in[i]*in[i])/4.0;
		if (out[i] < 0)
			out[i] = 0;
		out[i] = sqrt(out[i]);
	}
}

void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);

#if USE_GAUSSIAN_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian0,
		},
	};
#elif USE_TF_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
		},
	};
#else
	struct nlse_guess* guesses = NULL;
#endif


	//u32 bs[] = {4,8,12,16,24,32,48,64};
	//u32 os[] = {5,10,50,100,200,300,400,500,600,700,800,900,1000};
	//i64 os[] = {ORDER, 2*ORDER, 3*ORDER, 4*ORDER , 5*ORDER, 6*ORDER, 7*ORDER, 8*ORDER, 9*ORDER};
	i64 os[] = {1*ORDER, 2*ORDER, 3*ORDER, 4*ORDER, 5*ORDER};
	f64 g0s[] = {-1.0/ORDER, 1.0/ORDER};
	//i64 os[] = {9*ORDER};
	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation_gauss,
		.max_iterations = 1000,
		.max_integration_evals = 1e5,
		.error_tol = 1e-10,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.hamiltonian_mixing = 0.4,

		.zero_threshold = 1e-10,
		.gk=gk15
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(g0s)/sizeof(g0s[0]); ++j) {
		g0[0] = g0s[j];

		for (u32 i = 0; i < sizeof(os)/sizeof(os[0]); ++i) {
			u32 o;
			if (g0s[j] < 0)
				o = os[sizeof(os)/sizeof(os[0]) - i-1];
			else
				o = os[i];
			occupations[0] = o;

			sbmf_init();

			/* GP */
			struct nlse_result gp_res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
			f64 Egp = full_energy(settings, gp_res.coeff_count, component_count, gp_res.coeff, occupations, g0);

			/* BestMF */
			struct bestmf_result bmf_res = best_meanfield(settings, occupations[0], g0[0], guesses);

			/* RSPT */
			struct pt_result rs_res = rayleigh_schroedinger_pt_rf(settings, gp_res, 0, g0, occupations);

			/* ENPT */
			struct pt_result en_res = en_pt_rf(settings, gp_res, 0, g0, occupations);

			{
				FILE* fd = fopen("out", "a");
				fprintf(fd, "%lf\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
						g0[0]*(o-1),
						(Egp),
						(bmf_res.energy) - Egp,
						(rs_res.E0 + rs_res.E1 + rs_res.E2) - Egp,
						(rs_res.E0 + rs_res.E1 + rs_res.E2 + rs_res.E3) - Egp,
						(en_res.E0 + en_res.E1 + en_res.E2) - Egp,
						(en_res.E0 + en_res.E1 + en_res.E2 + en_res.E3) - Egp
						);
				fclose(fd);
			}

			sbmf_shutdown();
		}
	}

}
