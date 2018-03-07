#pragma once
#include "Node.h"
#include "Element.h"
#include "Local_Element.h"
class Grid {

public:
	double A;		//wysokoœæ siatki
	double B;		//szerokoœæ siatki

	int nA;		//liczba wêz³ów na œcianie A
	int nB;		//liczba wêz³ów na œcianie B

	int num_nodes;
	int num_elem;

	//int p; //liczba punktów ca³kowania

	Node* tab_nodes;
	Element* tab_elem;

	Local_Element* lE;
	double t0;
	double t_otoczenia;
	double max_czas;
	double krok_czasowy;
	double alfa;
	double cp; 
	double k;
	double ro;
	double* Det_Jakobian;   //tablica wyznaczników jakobianu - dla kazdego punktu calkowania osobno zapisujemy wyznacznik

	Grid();
	Grid(double, double, int, int, int, int, double, double, double, double, double, double, double, double);

	Grid(int, int, int, int, int, int, Node*, Element*, Local_Element, double, double,	double, double, double, double, double, double, double*);
	~Grid();

	void Ustaw_LH();		//liczy lokalne macierze[H] dla wszystkich elementów i wywo³uje funkcje set_LH(LH) dla ka¿dego elementu
	void Ustaw_LP();

	//void Utworz_wezly(Node*);
	void Utworz_wezly();
	//void Utworz_elementy(Element*);
	void Utworz_elementy();

	double** Oblicz_Jakobian(int id_elem, int p);

	double DetJ;

	double** Oblicz_LH(int id_elem);
	double* Oblicz_LP(int id_elem);
	double** Oblicz_C();

	double** oblicz_macierze(double** A, double** B);		//mno¿enie macierzy
	double** wektorXwektor(double* A, double* B);		//mno¿enie macierz x macierzy transponowanej
	double** MxM(double** A, double** B);







};