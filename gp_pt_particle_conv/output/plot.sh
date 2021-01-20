#!/bin/sh

	#set term epslatex;
	#set output 'imgs/gp_pt_particle_conv.eps';
gnuplot -p -e "
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set format y '%.3e';

	set xlabel '\$\\lambda\$';
	set ylabel 'Energy (a.\,u.)';
	plot
		'out' u 1:3 w linespoints lw 2 title 'Ebmf',
		'out' u 1:4 w linespoints lw 2 title 'rs2',
		'out' u 1:5 w linespoints lw 2 title 'rs3',
		'out' u 1:6 w linespoints lw 2 title 'en2',
		'out' u 1:7 w linespoints lw 2 title 'en3';
"
