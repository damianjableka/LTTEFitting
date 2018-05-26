#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define FALSE 0
#define TRUE 1
 //21.02.2012 12:17:13 dodano odczyt pliku z UP
 //21.02.2012 14:26:49 dziala poprawnie fcja Y i epo przetestowane dla znanych elementow ale z rysunku i z teorii sie zgadza
 //21.02.2012 14:36:50 prawidlowo odczytuje plik z wiakszo iloscia parametrow i ich granice
 //21.02.2012 15:08:09  zmieniona funkcja na ta do LITE i OC dodane rysowanie calych rozkladow i wywrsowanie wykresu z najmniejszym chikw
 //21.02.2012 15:46:22 dodano zapisywanie parametrow granicy oraz zamkniece wszystkich otwartych plikow
 //23.02.2012 09:44:01 zmienilem sposob uwzgledniania wagi w chikw bylo zle bo tylko jeden element byl dzielony przez nia, pomnozylem przez ta wage oba
 //23.02.2012 09:44:53 zmienilem opis pliku wejsciowego
 //23.02.2012 09:53:36 po zmianie z 9:44:01 program dopasowal idelana fcje, zmienilem sposob wyswietlania zawerzanych granic i same zawerzanie poszerzylem o +-10% przedzialu
 //23.02.2012 10:06:23 to zawerzenia o +-10% przedzialu powoduje ze graanice sie poszerzaja zamiast zawerzac wiec zobacyzmy co bez tego
 //23.02.2012 10:56:33 dodalem wyrysowywanie krzywej modelowej LITEqE^2 i samego qE^2ze skokiem co 10 epok 
 //24.02.2012 10:24:38 dodalem jeszce jeden plik na wejscie w nim musza sie znalesc ilosci losowan i powtorzen
 //24.02.2012 13:42:00 dodalem pelna fcje kwadratowa w postacji kanonicznej
 //28.02.2012 10:21:44 dodane inicjowanie ziarna generatora liczb losowych przed kazda zmiana grnaic
 //02.03.2012 09:03:06 dodalem obliczenie czasu do konca powtarzania granic
 //02.03.2012 13:20:10 dodalem wypisywanie residuow, zredukowanych danych na parabole, modelu bez praboli
 //02.03.2012 15:43:47 zmienilem format podawania q jako liczyby bez e-10 jest to zaszyte w modelu
 //02.03.2012 15:54:16 zmienilem z for po granicach na na while  chikw_min<0.3
 //02.03.2012 16:20:17 wywalilem co niepotrzebne wstawilem troche kometarza
 //05.03.2012 10:14:30 zmienile w pliku mcpar granice na zakres chikw
 //05.03.2012 10:14:51 usunaolem zapis do plikow rozkladow srand i granic bo to spowalnia raczej
 //20.03.2012 13:12:06 dodalem wypisanie na poczatku danych aby zobaczyc czy dobrze czyta
 //22.03.2012 12:17:16 zmienilem warunek konca na roznica chikw min i chikw max < Granica a możnaby jeszcze  zrobic chikw max > chikwmin(1+granica)
 //02.04.2012 13:17:42 dodalem rysowanie co jakis czas aktualnie co 30 zawerzen granic
 //02.04.2012 14:13:19 dodalem kolejny parametr w mcpar co ile ma drukowac wynik
 //03.04.2012 15:18:49 zmianilem ilosc przyjmowanych parametrow na >=4
 //03.04.2012 15:51:48 dodalem zapis calej bazy losowan, i jej odczyt jesli jest 5 argument jesli nie to pocztkowe losowania sa generowane
 //03.04.2012 15:52:25 dodalem zapis kopi parametrow poczatkowaych i nadpisywanie parametrow poczatkowych co iles zwerzen
 //04.04.2012 11:24:13 od teraz dziala zatrzymanie w dowolnym momencie i wznowienie symulacji dalej z nowymi parametrami poza iloscia wylosowanych pozycji w rozkladzie
 //13.04.2012 13:20:18 dodalem wyrysowywanie wykresu po kazdeych iteracjach przez system("gnuplot \'multiplot\'")
 //15.11.2012 11:26:17 plik wejsciowy sfromatowany jak w symulacji w matlabie zaiast z UP moze by zrpbic wybr najpierw
 //19.11.2012 16:02:19 kolejna proba zmianyz postaci kanonicznej na zwyklom aby nie dublowac parametro to si ei tak powinno ostatecznie dopasowac po zawerzeniu
 //10.12.2012 15:21:01 dodalem wybor w czwartym parametrze z rodejme pliku jaki podano 1 UP 0 matlab
 //06.02.2015 10:18:40 bootstrap a raczej jackknife z jednym elementem, 
// todo wszystko zapisuje sie do pliku jack na koncu po jego zamknieciu nalezy jeszcze policzyc srednia i stddev kazdej z uzyskanych wartosci.
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

struct zmienne //struktura jej zaleta jest 
{
 char nazwa[50];                // pole struktury odwolanie do pola p1.x=
 double dol; 			//limit dolny
 double gora;             	//limit gorny
};     

typedef struct zmienne Zm;   //definiuje nowy typ w c definiujemy zmienna typu particle: Particle p1, nie terzba tego robic ale wtedy usielibysmy defionwc struct s_particle p1

struct data //struktura zawierajaca dane 
{
 double jd;                // pole struktury odwolanie do pola p1.x=
 double ps; 			// odczytany z pliku pri/sec 1/0
 double w;             	//dczytana pliku z UP waga
};     

typedef struct data Da;   //definiuje nowy typ w c definiujemy zmienna typu particle


// na wejscie nazwa pliku z danymi i z zakresami zmiennych


void RandomInitialise(int ij,int kl) //plaski genberator liczb losowych inicjalizacja z biblioteki randomlib.c
{
   double s,t;
   int ii,i,j,k,l,jj,m;
   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
		ij = 1802;
		kl = 9373;
   }
   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);
   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }
   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = TRUE;
}

/* 
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
double RandomUniform(void)  
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test) 
   	RandomInitialise(1802,9373);
   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;
   return(uni);
}

double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}

int cmp( const void*a, const void*b)
{

 double one = *(double*) a;
   double two = *(double*)b;
if(one < two) return -1;
if(one > two) return +1;
return 0;
 }

double F(double x,double prawa, double e)            //funkcja opisujaca rownaie keplera
 {
  return(x-e*sin(x)-prawa);
 }

double epoka(double JD, double JD0, double P, double ps) //funkcja obliczajaca epoke dla danego JD i parametrow linowych
{
double epo=floor((JD-JD0)/P);
epo=((ps==1)?epo:(epo+0.5));  //z plikow z UP wyciagamy informacje o tym cyz min bylo pri (ps=1) czy sec (ps=0) podczas odczytu pliku z danymi
epo=((fabs(JD-(JD0+P*epo))>P/2)?epo+1:epo);
return epo;
}

double Y(double JD, double JD0, double P, double E)//fcja liczaca wartosc O-C z liniowych efemeryd
{
return(JD-(JD0+P*E));
}

double qLITE(double JD,double JD0,double P,double q,double p3,double A,double ecc,double T0,double w,double p,double r) //fcja opisujaca LITE+qE^2
{
double x_min=0; //dolna granica calokowania rownania keplera
double x_max=2.0* M_PI; //gorna granica calokowania rownania keplera
double tmp;
//tutaj mopze wyjsc pewnien problem z tymi epokami a JD trzeba obliczyc epoke i z niej znow jd
double ps=1;

double Epo=epoka(JD,JD0,P,ps);

JD=JD0+P*Epo;

double prawa=modf((JD-T0)/p3,&tmp)*M_PI*2; // prawa strona równania keplera
prawa=(prawa<0) ? prawa+2*M_PI : prawa;
double E;                                 
//zmienijszylem z 000001 do 0001 dla sprawdzenia czy to nie powoduje zwiechy
 while((x_max-x_min)>0.0001)          // metoda bisekji szukania miejsca zerowego rowniania keplera
  {
      E=(x_max+x_min)/2;
      x_max=(F(x_min,prawa,ecc)*F(E,prawa,ecc)<0) ? E : x_max;
      x_min=(F(x_max,prawa,ecc)*F(E,prawa,ecc)<0) ? E : x_min; 
  }

double nu =2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(E/2.0));
double DT=A/sqrt(1.0-pow(ecc*cos(w*M_PI/180.0),2.0))*((1.0-pow(ecc,2))/(1+ecc*cos(nu))*sin(nu+w*M_PI/180.0)+ecc*sin(w*M_PI/180.0));

// return(q*pow((Epo-p),2)*1e-10+r+DT); 19.11.2012 16:03:54 
return(q*pow((Epo),2)*1e-10+DT);
}

int main(int ile, char *tab[])
{
   FILE *dane;
   FILE *parametry;
   FILE *rozklad;
   FILE *wykres;
   FILE *los;
   FILE *model;
   FILE *par_mc;
   FILE *parametry_old;
   FILE *jack_f;
   FILE *wykres_j;
   FILE *model_j;
   FILE *out;
   double tmp,min,max;
   char ttmp[108];
   char pri[3];
   int c,d,e,ile_par,ile_danych,f,i,Ile_los,Ile,traf,Coile; //,t,dt;
   double chikw,Granice; //,czas;
int jaki_plik=0;
Zm par[11]; //tworze tablice zmienncyh typu Zm 
//---------------poczatek odczytu plikow-----------------------


 if(ile>=5)
  {
   if(fopen(tab[3],"r"))
   {
    par_mc=fopen(tab[3],"r");
    printf("odczytuje plik z danymi %s\n",tab[3]); 
    fscanf(par_mc,"%d %d %lf %d",&Ile_los,&Ile,&Granice,&Coile);  //trzeba zamienic Granice z int na double i w pliku wpisywac limit chikw zamiast ilosci i to wykorzystac na koncu

   if(fopen(tab[1],"r"))
   {
    dane=fopen(tab[1],"r");
    printf("odczytuje plik z danymi %s\n",tab[1]); 
  
    d=0;

      sscanf(tab[4],"%d",&jaki_plik); //26.04.2013 13:10:28 

     if(jaki_plik==1)// dodano 2012-12-10 15:20:28 wybor w 4 prametrze jaki rodzaj pliku podano 
     {
      printf("czytam plik z UP\n");
      while(fscanf(dane,"%108c\n",ttmp)!=EOF) // licze dlugosc pliku z danymi
      d++; 
     }
     else
     {
     printf("czytam plik prosty\n");
     while(fscanf(dane,"%lf%lf%3s\n",&tmp,&tmp,ttmp)!=EOF) //to dla pliku jak z matlaba 15.11.2012 11:46:42 
     d++;
     } 
    printf("plik ma %d linii\n",d);
   rewind(dane); // cofa sie na pocztake pliku
int ilosc_danych=d;
  Da dane_t[(d)]; //definieuje tablice danych typu Da o dlugosci d   

//------------ poczatek bootstrapa albo jackknife'a
printf("iolosc danyc : %d\n",ilosc_danych);
int jack;
double jack_a[ilosc_danych][11];
double dop_j[ilosc_danych][2];
int dpt;
 for(dpt=0;dpt<11;dpt++)
  {
  dop_j[dpt][0]=0;
  dop_j[dpt][1]=0;
  }

    d=ilosc_danych;
 
    for(e=0;e<d;e++)   //wczytuje plik z danymi do tablicy dane_t
     {
      if(jaki_plik)// dodano 2012-12-10 15:20:28 wybor w 4 prametrze jaki rodzaj pliku podano  zm 26.04.2013 13:10:17 
      {
       fscanf(dane,"%lf%12lf%8c%3s%9c%lf%68c\n",&tmp,&dane_t[e].jd,ttmp,pri,ttmp,&dane_t[e].w,ttmp); //plik z UP tak sformatowane na 108 kolumn
      }
      else
      {
      fscanf(dane,"%lf%lf%3s\n",&dane_t[e].jd,&dane_t[e].w,pri); //plik jak w symulacjach w matlabie z dodanymi pri albo sec w ostniej kolumnie 15.11.2012 11:33:15 
      dane_t[e].jd+=2400000.0; //to poniewarz w pliku z up jest 7 cyfr jd a w matlabie 5
      }
      dane_t[e].ps=(!strcmp(pri,"pri")?1:0);
     // printf("%lf\t%lf\t%lf\n",dane_t[e].jd,dane_t[e].ps,dane_t[e].w);
     }
printf("--wszystkie--\n");
      ile_danych=e;
int b;
int k=0;

ile_danych=e;

    if(k)
     printf("odczytano! %d wierszy\n",k);


    if(fopen(tab[2],"r")) //oczytuje plik z wartosciami i granicami
    {



    char nazwa[50]; //stworzenie nazwy do pliku kopii parametrow wejsciowych
   strcpy(nazwa,tab[2]);
     strcat(nazwa,".old");

  
    parametry=fopen(tab[2],"r");
  

//    d=0;
//    while(fscanf(parametry,"%s %lf %lf\n",ttmp,&tmp,&tmp)!=EOF) // licze dlugosc pliku z paramtrami = ilosc zmiennych
 //    d++;

//    rewind(parametry); //cofa wskaznik na oczatek pliku

//    Zm par[d]; //tworze tablice zmienncyh typu Zm 

     

    for(e=0;e<10;e++)   //wczytuje plik z szukanymi i ich granicami do tablicy par
     {
     fscanf(parametry,"%s %lf %lf\n",par[e].nazwa,&par[e].dol,&par[e].gora);
    }
    ile_par=e; //wpisujemy ilosc parametrow doasowania do zmiennej ile_par
    if(e)
     printf("odczytano i zapisano kopie do pliku %s! %d \n",nazwa, e);
fclose(parametry);


for(e=0;e<ile_par;e++)//odczyt parametrow
 printf("%s\t%lf\t%lf\n",par[e].nazwa,par[e].dol,par[e].gora);



double ep;


//jack dopisanie kolejnej lini z gory rozkladu do pliku jack
jack_f=fopen("jack","r");
rewind(jack_f);
int uu;
for( uu=0;uu<ilosc_danych;uu++)
  {
fscanf(jack_f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&jack_a[uu][1],&jack_a[uu][2],&jack_a[uu][3],&jack_a[uu][4],&jack_a[uu][5],&jack_a[uu][6],&jack_a[uu][7],&jack_a[uu][8],&jack_a[uu][9],&jack_a[uu][10],&jack_a[uu][0]);
printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",jack_a[uu][1],jack_a[uu][2],jack_a[uu][3],jack_a[uu][4],jack_a[uu][5],jack_a[uu][6],jack_a[uu][7],jack_a[uu][8],jack_a[uu][9],jack_a[uu][10],jack_a[uu][0]);
printf("%d\n",uu);
   for(dpt=0;dpt<11;dpt++)
   dop_j[dpt][0]+=jack_a[uu][dpt]/((double)ilosc_danych);
 }
fclose(jack_f);

    }
    else
    {
     printf("nie moge otworzyc pliku z parametrami, prawdopodobnie zla nazwa\n");
    }

   
// tu cala obsluga konca jacka

int uu=0;
out=fopen("mmc.out","w");
for(dpt=0;dpt<11;dpt++)
{
 for(uu=0;uu<ilosc_danych;uu++)
 {
  dop_j[dpt][1]+=pow((dop_j[dpt][0]-jack_a[uu][dpt]),2);
 }
dop_j[dpt][1]=sqrt(dop_j[dpt][1]/ilosc_danych);
fprintf(out,"%s \t %lf\t +-\t %lf \t(%.2lf %%)\n",((dpt>0)?par[dpt-1].nazwa:"chi^2"),dop_j[dpt][0],dop_j[dpt][1],fabs(dop_j[dpt][1]/dop_j[dpt][0]*100.));
printf("%s \t %lf\t +-\t %lf \t(%.2lf %%)\n",((dpt>0)?par[dpt-1].nazwa:"chi^2"),dop_j[dpt][0],dop_j[dpt][1],fabs(dop_j[dpt][1]/dop_j[dpt][0]*100.));
}

//rysunek ze srednich ale czy to ma sens??
wykres_j=fopen("wykres_j","w");
  // prawidlowo rysuje O-C !!  zarowno dla teroi fcja LITE jak i fcja Y dla danych
double ql, ep, oc;
	for(e=0;e<ile_danych;e++)
	 {
	  ql=qLITE(dane_t[e].jd,dop_j[1][0],dop_j[2][0],dop_j[3][0],dop_j[4][0],dop_j[6][0],dop_j[6][0],dop_j[7][0],dop_j[8][0],dop_j[9][0],dop_j[10][0]);
	  ep=epoka(dane_t[e].jd,dop_j[1][0],dop_j[2][0],dane_t[e].ps);
	  oc=Y(dane_t[e].jd,dop_j[1][0],dop_j[2][0],ep);
  	  fprintf(wykres_j,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",ep,ql,oc,dane_t[e].jd,dane_t[e].w,oc-ql,oc-(dop_j[0][3]*pow((double)ep,2)*1e-10));
         }
          printf("stworzylem plik \"wykres\" z danymi:\n epoka \t model \t O-C \t jd\t waga\t dane-model\t dane-q,p,r\n");

model_j=fopen("model_j","w");
double jd;
for(e=epoka(dane_t[0].jd,dop_j[1][0],dop_j[2][0],dane_t[0].ps);e<epoka(dane_t[ile_danych-1].jd,dop_j[1][0],dop_j[2][0],dane_t[ile_danych-1].ps);e+=10)
 {
  jd=dop_j[1][0]+dop_j[2][0]*(double)e;
  ql=qLITE(jd,dop_j[1][0],dop_j[2][0],dop_j[3][0],dop_j[4][0],dop_j[5][0],dop_j[6][0],dop_j[7][0],dop_j[8][0],dop_j[9][0],dop_j[10][0]);
  fprintf(model_j,"%lf\t%lf\t%lf\t%lf\n",(double)e,ql,(dop_j[3][0]*pow(((double)e),2)*1e-10),ql-(dop_j[3][0]*pow((double)e,2)*1e-10));
 }

fclose(wykres_j);
fclose(model_j);

system("gnuplot \'multiplot_j\'");



  fclose(dane);





   }
   else
   {
    printf("nie moge odczytac pliku z danymi, prawdopodobnie zla nazwa\n");
   }
   }
   else
   {
    printf("nie moge odczytac pliku z parametrami losowan, prawdopodobnie zla nazwa\n");
   }
  }
  else
 {
  printf("na wejsciu podaj nazwe pliku z  UP z danymi\nn\tjd\t8xchar\tpri/sec\t9xchar\twaga\t 68xchar\n\n a nastepnie nazwe pliku z parametrami i zakresami w formacie:\n nazwa_zmiennej\t dony_limit\t gorny_limit\n\n oraz plik z parametrami losowania montecarlo w formacie:\n ilosc_losowan_wstepnych \t ilosc_losowan_kolejnych\t warunek konca delta ochylenia \t ilosc_zawerzen_granic\n\n rodzaj pliku \n1-UP \n0-prosty (jak w matlabie z pri/sec w 3 kolumnie)\n \n opcjonalnie nazw pliku z rozkladem jesli chcesz kontynuwoac obliczenia\n");
 }
return(0);
}

