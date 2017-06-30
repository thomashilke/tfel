set term epslatex size 11.25cm,6.95cm color colortext 10
set output 'quadrature-convergence.tex'

set xlabel '$h$'
set ylabel 'Erreur $\displaystyle\frac{\abs{I - L_h(f)}}{\abs{I}}$'

set xtics format '%3.1e'
set ytics format '%3.1e'

set key top left
set key spacing 1.25

set logscale xy

plot '../data/quadrature_1d_convergence.dat' u 1:2 w lp title 'Gauss 1 point', \
     '' u 1:3 w lp title 'Gauss 2 point', \
     '' u 1:4 w lp title 'Gauss 3 point', \
     '' u 1:5 w lp title 'Gauss 4 point', \
     '' u 1:6 w lp title 'Gauss 5 point'
