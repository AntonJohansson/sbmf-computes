#include <sbmf/sbmf.h>

#include <stdio.h>

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define USE_RANDOM_GUESS 1

#define PLOT 0

#if PLOT
	#include <plot/plot.h>
#endif













void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
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
#elif USE_RANDOM_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
#else
	struct nlse_guess* guesses = NULL;
#endif

	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.error_tol = 1e-14,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.orbital_mixing = 0.0,
		.hamiltonian_mixing = 0.5,

		.orbital_choice = NLSE_ORBITAL_LOWEST_ENERGY,

		.gk=gk20
    };

	OMEGA = 0.01;

	const u32 component_count = 2;

	i64 Ns[] = {4, 8, 16, 32, 64, 128, 256, 512};

	FILE* fd = fopen("out", "a");
	fprintf(fd, "# N\tE\tRS2\tRS3\tEN2\tEN3\n");

	for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {
		i64 N = Ns[i];
		i64 occupations[] = {N,N};
		f64 g0[] = {
			 1.0/((f64)N-1), -0.5/((f64)N-1),
			-0.5/((f64)N-1),  1.0/((f64)N-1)
		};

		sbmf_init();
		struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
		f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
		printf("\nfull energy: %lf\n", Efull);

		struct pt_result rspt = rayleigh_schroedinger_pt_rf_2comp(settings, res, g0, occupations);
		struct pt_result enpt = en_pt_2comp(settings, res, g0, occupations);
		fprintf(fd, "%ld\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
			N,
			Efull,
			rspt.E0+rspt.E1+rspt.E2,
			rspt.E0+rspt.E1+rspt.E2+rspt.E3,
			enpt.E0+enpt.E1+enpt.E2,
			enpt.E0+enpt.E1+enpt.E2+enpt.E3
			);
		sbmf_shutdown();
	}

	fclose(fd);

#if PLOT
    {
        const u32 N = 256;
        plot_init(800, 600, "gp2c");
        f32 potdata[N], adata[N], bdata[N];
        sample_space sp = make_linspace(1, -5, 5, N);

        for (u32 i = 0; i < N; ++i) {
            f64 x = sp.points[i];
            potdata[i] = (f32) ho_potential(&x,1,0);
        }
        push_line_plot(&(plot_push_desc){
                .space = &sp,
                .data = potdata,
                .label = "potential",
                });


        f64 sample_in[N];
        for (u32 i = 0; i < N; ++i) {
            sample_in[i] = (f64) sp.points[i];
        }

        for (u32 i = 0; i < res.component_count; ++i) {
            f64 sample_out[N];
            ho_sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

            f32 data[N];
            for (u32 k = 0; k < N; ++k) {
                data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
            }
            push_line_plot(&(plot_push_desc){
                    .space = &sp,
                    .data = data,
                    .label = plot_snprintf("comp: %u", i),
                    .offset = res.energy[i],
                    });
        }
		plot_update_until_closed();
        plot_shutdown();
	}
#endif


}
