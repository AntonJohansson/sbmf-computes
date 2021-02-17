#include <sbmf/sbmf.h>

#include <stdio.h>

#define NA 4
#define NB 4

//#define GAA (1.0/1e5)
//#define GAA (-4.0)
//#define GAA (-10.0/(NA-1))
//#define GAA (0.5/(NA-1))
#define GAA (1.0/3.0)
#define GAB (-1.0/6.0)
#define GBA (-1.0/6.0)
#define GBB (1.0/3.0)
//#define GAB (-1.0/(NB))
//#define GBA (-1.0/(NA))
//#define GBB (0.5/(NB-1))

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define COMPONENT_COUNT 2

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

void tf(f64* out, f64* in, u32 len, void* data) {
	f64 mu = 0.5 * sqrt(3*GAA*(NA-1)/2);
	for (u32 i = 0; i < len; ++i) {
		out[i] = (mu - 0.5*in[i]*in[i])/4.0;
		if (out[i] < 0)
			out[i] = 0;
		out[i] = sqrt(out[i]);
	}
}

void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i) {
		//out[i] = (1 - 0.5 * in[i]*in[i]);
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);// + gaussian(in[i] - 1.0, 0.0, 0.2);
	}
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}


void expnx(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = exp(-fabs(in[i]));
	}
}


int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();

	OMEGA = 0.1;

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
#elif USE_TF_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
		},
		[1] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
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

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.error_tol = 1e-14,

		.num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.hamiltonian_mixing = 0.5,

		.orbital_choice = NLSE_ORBITAL_LOWEST_ENERGY,

		.gk=gk20
    };

	const u32 component_count = COMPONENT_COUNT;

	struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
	f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
	printf("\nfull energy: %lf\n", Efull);

	nlse_write_to_binary_file("outbin", res);

#if 1
	{
		//struct pt_result ptres = rayleigh_schroedinger_pt_rf(settings, res, 0, g0, occupations);
		//struct pt_result ptres = rspt_1comp_cuda_new(&settings, res, 0, g0[0], occupations[0]);
		//struct pt_result ptres = enpt_1comp_cuda_new(&settings, res, 0, g0[0], occupations[0]);
		//struct pt_result ptres = enpt_1comp_cuda(settings, res, 0, g0, occupations);
		struct pt_result ptres = rspt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
		//struct pt_result ptres = enpt_2comp_cuda_new(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
		//struct pt_result ptres = en_pt_2comp(settings, res, g0, occupations);
		printf("E0:          %.15lf\n", ptres.E0);
		printf("E1:          %.15lf\n", ptres.E1);
		printf("E2:          %.15lf\n", ptres.E2);
		printf("E3:          %.15lf\n", ptres.E3);
		printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);
		printf("diff: %.15lf\n", Efull - (ptres.E0+ptres.E1+ptres.E2+ptres.E3));
	}
#endif

	sbmf_shutdown();
}
