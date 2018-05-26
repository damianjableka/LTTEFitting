# LTTEFitting
Monte Carlo fitting procedure to fit 8 parameter space to O-C

Tested only on ScientificLinux 5 - 7 
compilation:
gcc -Wall -lm mmc36beta.c -o mmc36beta

required c libraries:
stdio.h
string.h
math.h
stdlib.h
time.h

required Gnuplot > 4.0

par - Monte Carlo searching parameters boundaries
mcpar - Monte Carlo procedure parameters
run mmc* without any parameters to detailed discription
