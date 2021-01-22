#!/bin/sh

	#set term epslatex;
	#set output 'imgs/gp_pt_1c_VH_particle_conv.eps';
	#set output 'imgs/gp_pt_1c_VHG_particle_conv.eps';
gnuplot -p -e "
	set size 0.75,0.75;
	set key at graph 0.75,0.6;
	set format y '%.1f';
	set key spacing 1.5;

	set xlabel '\$\\lambda = 0.1\\cdot(N-1)\$';
	set ylabel 'Energy relative to \$E_\\mathrm{GP}\$ (a.\,u.)';

	plot
		'out_N_2c_VH' u 1:3 w linespoints lw 2 lc 4 dt 2 	title '\$E_\\mathrm{RSPT2}\$',
		'out_N_2c_VH' u 1:4 w linespoints lw 2 lc 4 		title '\$E_\\mathrm{RSPT3}\$',
		'out_N_2c_VH' u 1:5 w linespoints lw 2 lc 6 dt 2 	title '\$E_\\mathrm{ENPT2}\$',
		'out_N_2c_VH' u 1:6 w linespoints lw 2 lc 6 		title '\$E_\\mathrm{ENPT3}\$';
"