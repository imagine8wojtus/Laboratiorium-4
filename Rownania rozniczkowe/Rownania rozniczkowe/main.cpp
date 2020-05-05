#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rk4.h"

double lambda;                                  //jako zmienna globalna
double funkcja(double t, double y);

void main(){

	double h, tk, y0, t0;
	int N;

	printf("Podaj wartosc wspolczynnika lambda:\n");           //pobieranie poszczegolnych wartosci od uzytkownika
	scanf("%lf", &lambda);
	printf("Podaj warunek poczatkowy t0:\n");
	scanf("%lf", &t0);
	printf("Podaj warunek poczatkowy y0:\n");
	scanf("%lf", &y0);
	printf("Podaj wartosc tk:\n");
	scanf("%lf", &tk);

	FILE* f;
	f = fopen("rownania_rozniczkowe.csv", "w");   //wpisywanie do pliku 

	double blad_e, blad_rk;

	for (int k = 0; k <= 6; k++)                //petla dla wartosci N
	{
		N = pow(2.,(double)k);
		h = (tk - t0) / N;                     //z rownania (tk - t0)/N wyjdzie nam dlugosc jednego kroku

		double* t,* yan,* ye,* yrk;                     
		t = (double*)malloc((N+1) * sizeof(double));    //tworzenie tablic jednowymiarowych dla poszczegolnych wartosci, N+1 poniewaz mamy liczbe krokow plus dla t0
		yan = (double*)malloc((N+1)* sizeof(double));   //tworzymy je aby przechowywac w nich wyliczone wartosci
		ye = (double*)malloc((N+1) * sizeof(double));
		yrk = (double*)malloc((N+1) * sizeof(double));

		for (int i = 0; i <= N; i++) {
			t[i] = t0 + (h * i);
			yan[i] = y0 * exp((t[i] - t0) * lambda);       //wartosci wyliczane analitycznie       
		}

		ye[0] = y0;
		for (int i = 0; i < N; i++) {
			ye[i + 1] = ye[i] + h * funkcja(t[i], ye[i]);      //wartosci wyliczane metoda Eulera 
		}

		yrk[0] = y0;
		for (int i = 0; i < N; i++) {
			yrk[i + 1] = rk4(t[i], yrk[i], h, funkcja);        //wyliczanie wartosci metoda Rungego-Kutty poprzez wywolanie funckji rk4
		}

		blad_e = fabs((ye[N] - yan[N]) / yan[N]);               //wyliczanie bledu dla metody Eulera
		blad_rk = fabs((yrk[N] - yan[N]) / yan[N]);             //wyliczanie bledu dla metody Rungego-Kutty

		fprintf(f, "%d\t%lf\t%lf\t%lf\n", N, h, blad_e, blad_rk);

		free(t);
		free(yan);
		free(ye);
		free(yrk);

	 }
	fclose(f);
};

double funkcja(double t, double y)
{
	return lambda * y;
}


