echo """
set terminal png size 1920,1080
set output 'plots/bicubic.png'

set multiplot layout 2,3 rowsfirst


set pm3d map

# similar to matlab color palette
set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',\
                      4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000',\
                      8 '#7f0000')


set xlabel 'x'
set ylabel 'y'
#
# ----- EVAL graph -----

set title 'eval' font ',25'
splot 'out/bicubic.dat' us 1:2:3

# ----- EVAL x graph -----

set title 'eval x' font ',25'
splot 'out/bicubic.dat' us 1:2:4

# ----- EVAL y graph -----

set title 'eval y' font ',25'
splot 'out/bicubic.dat' us 1:2:5

# ----- EVAL xx graph -----

set title 'eval xx' font ',25'
splot 'out/bicubic.dat' us 1:2:6

# ----- EVAL yy graph -----

set title 'eval yy' font ',25'
splot 'out/bicubic.dat' us 1:2:7

# ----- EVAL xy graph -----

set title 'eval xy' font ',25'
splot 'out/bicubic.dat' us 1:2:8

""" | gnuplot
