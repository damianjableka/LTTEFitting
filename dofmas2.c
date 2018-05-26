#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define FALSE 0
#define TRUE 1
#define dniwroku 365.256


//log zmian
//2013-04-26 20:31:22 dziala poprawnie liczy mase z obliczonej funkcji mas 
//
//--------------------
int main(int ile, char *tab[])
{
   FILE *parametry;
   double tmp;
  double p3, A, ecc, w, M1, M2, ain, i12, Pin;
double G=6.67e-8; //dyn cm^2/g^2
double solarM=1.988435e33; // g
double AU=1.49597871e13; //cm
double sekwdniu=24.*3600.; // sec/dobe
 if(ile>1)
  {
   if(fopen(tab[1],"r"))
   {

    parametry=fopen(tab[1],"r");

   fscanf(parametry,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp,&Pin,&tmp,&p3,&A,&ecc,&tmp,&w);
    printf("p3 %lf\nA %lf\ne %lf\nw %lf\n",p3,A,ecc,w);


if(ile==5)
{
sscanf(tab[2],"%lf",&M1);
sscanf(tab[3],"%lf",&M2);
sscanf(tab[4],"%lf",&i12);
   ain=pow((G*(M1+M2)*solarM*pow(Pin*sekwdniu,2)/(4.*pow(M_PI,2))),1/3.)/AU;
printf("M1 %lf\nM2 %lf\nain %lf\ni12 %lf\n",M1,M2,ain,i12);
}
else
{
printf("M1\nM2\nain \ni12\n");
}
  }else
   {
     printf("nie udalo sie otworzyc pliku \n\n ");
    }
  }
  else
 {
  printf("na wejsciu podaj lokalizacje pliku rozklady\n oraz ewentualnie M1 M2 iin\n ");
 }
return(0);
}

