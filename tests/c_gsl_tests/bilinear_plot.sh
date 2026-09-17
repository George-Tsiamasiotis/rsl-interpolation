echo """
set terminal png size 1920, 540
set output 'out/bilinear.png'

set multiplot layout 1,3 rowsfirst


set pm3d map

# similar to matlab color palette
set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',\
                      4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000',\
                      8 '#7f0000')


set xlabel 'x'
set ylabel 'y'

# ----- EVAL graph -----

set title 'eval' font ',25'
splot 'out/bilinear.dat' us 1:2:3

# ----- EVAL x graph -----

set title 'eval x' font ',25'
splot 'out/bilinear.dat' us 1:2:4

# ----- EVAL y graph -----

set title 'eval y' font ',25'
splot 'out/bilinear.dat' us 1:2:5
""" | gnuplot
