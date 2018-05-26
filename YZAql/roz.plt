reset
set term postscript landscape enhanced color lw 1 dashed 'Helvetica' 14
set autoscale
set xl 'szerokość przedzialu'
set yl 'maksymalna różnica ilosc punktow w przedziale'
set output 'grupowanie.eps'
plot 'zgrupowanie' u 1:2
reset
