set terminal pngcairo size 1000,1000 enhanced font 'Arial,12'
set output 'fractal_completo.png'

set title "Fractal de Bairstow"
set xlabel "Chute Inicial r0"
set ylabel "Chute Inicial s0"

set xrange [-5:5]
set yrange [-5:5]

set datafile separator ","

set palette rgbformulae 22,13,-31

set cbrange [0:20]
unset colorbox

plot 'fractal_bairstow.csv' using 1:2:(int($4) % 20) with points pt 5 ps 0.5 palette notitle
