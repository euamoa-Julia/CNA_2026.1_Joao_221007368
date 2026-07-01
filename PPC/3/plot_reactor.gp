# =====================================================================
# SEÇÃO 1: CAMPOS DE CORES (MAPAS TÉRMICOS PM3D)
# =====================================================================
set terminal pngcairo size 900,650 enhanced font 'Verdana,10'
set pm3d map explicit
unset pm3d
set border 31 lt 1 lc rgb '#444444' lw 1
set xrange [0:0.01]
set yrange [0:300]
set cblabel "Temperatura (°C)" font 'Verdana,10,Bold'
# Paleta estilo "magma" (matplotlib), do preto/roxo escuro ao amarelo pálido
set palette defined ( \
    0.000 '#000004', \
    0.111 '#1c1044', \
    0.222 '#4b0c6b', \
    0.333 '#781c6d', \
    0.444 '#a52c60', \
    0.556 '#cf4446', \
    0.667 '#ed6925', \
    0.778 '#fb9a06', \
    0.889 '#f7d03c', \
    1.000 '#fcfdbf' )
# 1.1 Mapa
set output 'campo.png'
set title "Campo de Cores - Perfil Térmico" font 'Verdana,12,Bold'
splot 'mapa.dat' using 1:2:3 notitle with pm3d
unset output
# =====================================================================
# SEÇÃO 2: EVOLUÇÃO TEMPORAL (GRÁFICOS DE LINHAS T vs t)
# =====================================================================
reset  
set terminal pngcairo size 900,600 enhanced font 'Verdana,10'
set grid xtics ytics lc rgb '#dddddd' lt 1 lw 1
set key left top box lt 1 lc rgb '#bbbbbb'
set xlabel "Tempo (s)" font 'Verdana,10,Bold'
set ylabel "Temperatura (°C)" font 'Verdana,10,Bold'
# Mesmas 7 cores extraídas do colormap "magma", do mais escuro (nó central) ao mais claro (borda)
set linetype 1 lc rgb "#000004" lw 2.5  
set linetype 2 lc rgb "#221150" lw 2.5  
set linetype 3 lc rgb "#5f187f" lw 2.5  
set linetype 4 lc rgb "#982d80" lw 2.5  
set linetype 5 lc rgb "#d3436e" lw 2.5  
set linetype 6 lc rgb "#f8765c" lw 2.5  
set linetype 7 lc rgb "#febb81" lw 2.5  
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
