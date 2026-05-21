# =====================================================================
# SEÇÃO 1: CAMPOS DE CORES (MAPAS TÉRMICOS PM3D)
# =====================================================================
set terminal pngcairo size 900,650 enhanced font 'Verdana,10'

set pm3d map explicit
set pm3d noborder
set border 31 lt 1 lc rgb '#444444' lw 1
set xrange [0:0.01]
set yrange [0:300]
set cblabel "Temperatura (°C)" font 'Verdana,10,Bold'

set palette defined ( 0.0 "#00008b", 0.2 "#3b4cc0", 0.4 "#8dd3c7", 0.6 "#ffffb3", 0.8 "#f7745b", 1.0 "#b40426" )

# 1.1 Mapa
set output 'campo.png'
set title "Campo de Cores - Perfil Térmico" font 'Verdana,12,Bold'
splot 'mapa.dat' using 1:2:3 notitle with pm3d

# =====================================================================
# SEÇÃO 2: EVOLUÇÃO TEMPORAL (GRÁFICOS DE LINHAS T vs t)
# =====================================================================
reset  

set terminal pngcairo size 900,600 enhanced font 'Verdana,10'
set grid xtics ytics lc rgb '#dddddd' lt 1 lw 1
set key left top box lt 1 lc rgb '#bbbbbb'

set xlabel "Tempo (s)" font 'Verdana,10,Bold'
set ylabel "Temperatura (°C)" font 'Verdana,10,Bold'

set linetype 1 lc rgb "#b40426" lw 2.5  
set linetype 2 lc rgb "#e7745b" lw 2.5  
set linetype 3 lc rgb "#f7b89c" lw 2.5  
set linetype 4 lc rgb "#dddcdc" lw 2.5  
set linetype 5 lc rgb "#aac7fd" lw 2.5  
set linetype 6 lc rgb "#6f92f3" lw 2.5  
set linetype 7 lc rgb "#3b4cc0" lw 2.5  

# 2.1 Curvas
set output 'curvas.png'
set title "Evolução Temporal da Temperatura" font 'Verdana,12,Bold'
plot 'linhas.dat' using 1:2 title 'Nó 1 (Centro)' with lines lt 1, \
     '' using 1:3 title 'Nó 2' with lines lt 2, \
     '' using 1:4 title 'Nó 3' with lines lt 3, \
     '' using 1:5 title 'Nó 4' with lines lt 4, \
     '' using 1:6 title 'Nó 5' with lines lt 5, \
     '' using 1:7 title 'Nó 6' with lines lt 6, \
     '' using 1:8 title 'Nó 7 (Borda)' with lines lt 7

unset output
