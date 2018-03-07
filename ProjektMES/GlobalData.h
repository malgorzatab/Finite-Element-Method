#pragma once
#include "Grid.h"
class GlobalData {

public:
	int num_nodes;		//16 nA*nB
	int num_elem;		//9  (nA-1)(nB-1)

	double** global_H; //zale�ne od ilo�ci w�z��w
	double* global_P;

	double* temperatury;

	GlobalData();
	GlobalData(int, int, double**, double*, double*);
	//Grid* grid;
	void wczytaj();
	double* gauss(double** A, double* B);
	void oblicz_GH(Grid grid);
	void oblicz_GP(Grid grid);

	void utworz_Siatke();

	double t0; 
	double max_czas; 
	double time_step; 
	double t_otoczenia;
	double alfa;
	double A; 
	double B; 
	int nA; 
	int nB;
	double cp; 
	double k; 
	double ro; 

};


//tworze grid, wpisuje z pliku, w gridzie tworz� si� w�z�y i elementy, ustawiam temp pocz�tkowe, p�tla po elementach 