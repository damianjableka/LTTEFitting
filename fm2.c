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
struct zmienne //struktura jej zaleta jest 
{
 char nazwa[50];                // pole struktury odwolanie do pola p1.x=
 double wartosc; 			//limit dolny

};     

typedef struct zmienne Zm; 


int main(int ile, char *tab[])
{
   FILE *parametry;
   double tmp;
   char ttmp[108];
  double p3, A, ecc, w, M1, M2, el;
double m3,u,n, a12, a3, ain, Q, Qst, i12, i312, st;
 if(ile>1)
  {
   if(fopen(tab[1],"r"))
   {

parametry=fopen(tab[1],"r");
    printf("odczytuje plik z parametrami %s\n",tab[1]); 

   fscanf(parametry,"%s %lf \n%s %lf \n%s %lf \n%s %lf \n%s %lf \n%s %lf \n%s %lf \n%s %lf\n",ttmp,&p3,ttmp,&A,ttmp,&ecc,ttmp,&w,ttmp,&M1,ttmp,&M2,ttmp,&ain,ttmp,&i12);
    printf("p3 %lf\nA %lf\necc %lf\nw %lf\nM1 %lf \nM2 %lf\nain %lf\ni12 %lf\n",p3,A,ecc,w,M1,M2,ain,i12);
     double f=1/pow(p3/dniwroku,2)*pow(((173.15*A)/sqrt(1-pow(ecc,2)*pow(cos(w*M_PI/180.),2))),3);
     printf("fcja mas = %lf\n",f);
//for(el=0.2;el<=1;el+=0.01)
printf("sin(i)\t\t i \t\t m3\t\t a12 \t\t a3 \t\t Q\t\t Qst\t\t i312\t\t st\n");
for(el=15.0;el<90;el+=1)
{
u=sin(el*M_PI/180); //sin i
n=M1+M2;

double wyraz=pow(2.*pow(f,3)+18.*pow(f,2)*n*pow(u,3)+3.*sqrt(3.)*sqrt(4*pow(f,3)*pow(n,3)*pow(u,9)+27.*pow(f,2)*pow(n,4)*pow(u,12))+27.*f*pow(n,2)*pow(u,6),1./3.);
m3=wyraz/(3*pow(2.,1./3.)*pow(u,3));
m3-=(pow(2.,1./3.)*(-pow(f,2)-6*f*n*pow(u,3)))/(3*pow(u,3)*wyraz);
m3+=f/3*pow(u,2);

a12=173.15*A/sqrt(1-pow(ecc,2)*pow(cos(w*M_PI/180),2))/u;
a3=a12*m3/n;
i312=(asin(u)*180/M_PI-i12);
Q=a3*(1-ecc)/ain;
Qst=3.*pow((1+m3/n),1./3.)*pow(1-ecc,-1./6.)*pow((7./4.+cos(i312*M_PI/180.)/2-pow(cos(i312*M_PI/180.),2)),1./3.);
st=1;
if(Q<Qst)
st=0;
printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t%lf\t%lf\n",u,asin(u)*180/M_PI,m3,a12,a3,Q,Qst,i312,st);
}
   }
   else
   {
     printf("nie udalo sie otworzyc pliku \n\n ");
    }
  }
  else
 {
  printf("na wejsciu podaj plik w formacie \np3\n A\n e\n w\n M1\n M2\n ain\n i12\n najpierw nazwa potem wartosc\n program dofmas generuje czesc tego pliku z pliku rozklad trzeba dodac do tak wygenerowanego jeszcze M1 M2 oraz ain i12");
 }
return(0);
}

