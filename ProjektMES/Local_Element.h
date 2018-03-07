#pragma once
class Local_Element {

public:
	double ksi[4];
	double eta[4];		//wspolrzedne punktow calkowania w objeto�ci w ukl lokalnym 

	double ksi_pow[8];
	double eta_pow[8];	//wspolrzedne punkt�w ca�kowania na powierzchni w ukl lokalnym


	double N_V[4][4];  //funkcje kszta�tu w punktach ca�kowania w obj�to�ci 4 punkty - 4 funkcje w ka�dym punkcie
	double N_P[8][4];  //funkcje kszta�tu w punktach ca�kowania po powierzchni 8 punkt�w, po 2 na powierzchnie- 4 funkcje dla ka�dego punktu

	double dNdKsi[4][4]; //4 punkty ca�kowania i dla ka�dego po 4 warto�ci pochodnych 
	double dNdEta[4][4];

	Local_Element();

	Local_Element(double* ksi, double*eta, double* ksi_pow, double* eta_pow, double**N_V, double** N_P, double** dNdKsi, double** dNdEta);
	~Local_Element();

	void setKsi();
	void setEta();

	void obliczN_V();
	
	void obliczN_P();
	
	void obliczdNdKsi();
	void obliczdNdta();

};