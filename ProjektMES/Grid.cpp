#include "Grid.h"
#include<iostream>
using namespace std;

Grid::Grid() {
	this->A = 0;
	this->B = 0;
	this->nA = 0;
	this->nB = 0;
	this->num_nodes = 0;
	this->num_elem = 0;
	this->tab_nodes = new Node[this->num_nodes];
	this->tab_elem = new Element[this->num_elem];
	this->lE = new Local_Element();
	this->t0 = 0.0;
	this->t_otoczenia = 0.0;
	this->max_czas = 0.0;
	this->krok_czasowy = 0.0;
	this->alfa = 0.0;
	this->cp = 0.0;
	this->k = 0.0;
	this->ro = 0.0;
	this->Det_Jakobian = new double[4];		

	for (int i = 0; i < 4; i++) {		
			this->Det_Jakobian[i] = 0;		
	}

	//Utworz_elementy(this->tab_elem);
	//Utworz_wezly(this->tab_nodes);
	Utworz_elementy();
	Utworz_wezly();
}

Grid::Grid(int A, int B, int nA, int nB, int num_nodes, int num_elem, Node* tab_nodes, Element* tab_elem, Local_Element lE, double t0, double t_otoczenia,
	double max_czas, double krok_czasowy, double alfa, double cp, double k, double ro, double* Jakobian) {
	this->A = A;
	this->B = B;
	this->nA = nA;
	this->nB = nB;
	this->num_nodes = num_nodes;
	this->num_elem = num_elem;
	this->tab_nodes = new Node[num_nodes];
	this->tab_elem = new Element[num_elem];
	this->lE = new Local_Element();
	this->t0 = t0;
	this->t_otoczenia = t_otoczenia;
	this->max_czas = max_czas;
	this->krok_czasowy = krok_czasowy;
	this->alfa = alfa;
	this->cp = cp;
	this->k = k;
	this->ro = ro;
	this->Det_Jakobian = new double[4];

	for (int i = 0; i < 4; i++) {
			this->Det_Jakobian[i] = 0;		
	}

	Utworz_elementy();
	Utworz_wezly();
}


Grid::Grid(double A, double B, int nA, int nB, int num_nodes, int num_elem, double t0, double t_otoczenia,
	double max_czas, double krok_czasowy, double alfa, double cp, double k, double ro) {
	this->A = A;
	this->B = B;
	this->nA = nA;
	this->nB = nB;
	this->num_nodes = num_nodes;
	this->num_elem = num_elem;
	this->tab_nodes = new Node[num_nodes];
	this->tab_elem = new Element[num_elem];
	this->lE = new Local_Element();
	this->t0 = t0;
	this->t_otoczenia = t_otoczenia;
	this->max_czas = max_czas;
	this->krok_czasowy = krok_czasowy;
	this->alfa = alfa;
	this->cp = cp;
	this->k = k;
	this->ro = ro;
	this->Det_Jakobian = new double[4];

	for (int i = 0; i < 4; i++) {
		this->Det_Jakobian[i] = 0;
	}

	Utworz_elementy();
	Utworz_wezly();
}

Grid::~Grid() {}

void Grid::Ustaw_LH() {
	for (int i = 0; i <this->num_elem; i++) {
		this->tab_elem[i].setlH(Oblicz_LH(this->tab_elem[i].id_element ));	//dla ka¿dego elementu oblicza macierz H
	}
}

void Grid::Ustaw_LP() {
	for (int i = 0; i <this->num_elem; i++) {
		this->tab_elem[i].setlP(Oblicz_LP(this->tab_elem[i].id_element));	//dla ka¿dego elementu oblicza macierz P
	}
}

double** Grid::Oblicz_LH(int id_elem) {
	
	//macierz H
	double** H = new double*[4];
	for (int l = 0; l < 4; l++) {
		H[l] = new double[4];
	}

	//1 czêœæ macierzy H
	double** macierzH1 = new double*[4];
	for (int i = 0; i < 4; i++) {
		macierzH1[i] = new double[4];
	}

	//2 czêœæ macierzy H
	double** macierzH2 = new double*[4];
	for (int i = 0; i < 4; i++) {
		macierzH2[i] = new double[4];
	}

	//3 czêœæ macierzy H
	double** macierzH3 = new double*[4];
	for (int i = 0; i < 4; i++) {
		macierzH3[i] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = 0;
			macierzH1[i][j] = 0;
			macierzH2[i][j] = 0;
			macierzH3[i][j] = 0;
		}
	}

	double** dNdXY = new double*[2];
	for (int k = 0; k < 2; k++) {
		dNdXY[k] = new double[4];
	}

	double* dNdX = new double[4];
	double* dNdY = new double[4];


	for (int p = 0; p < 4; p++) {				//petla po p ca³kowania

		dNdXY = Oblicz_Jakobian(id_elem, p);


		for (int i = 0; i < 4; i++) {
			dNdX[i] = 0;
			dNdY[i] = 0;

		}

		dNdX[0] = dNdXY[0][0];
		dNdX[1] = dNdXY[0][1];
		dNdX[2] = dNdXY[0][2];
		dNdX[3] = dNdXY[0][3];
		dNdY[0] = dNdXY[1][0];
		dNdY[1] = dNdXY[1][1];
		dNdY[2] = dNdXY[1][2];
		dNdY[3] = dNdXY[1][3];

		//transponowanie macierz dndX i dNdY, mno¿enie i sumowanie tych macierzy

		double** X = new double*[4];
		for (int i = 0; i < 4; i++) {
			X[i] = new double[4];
		}

		X = wektorXwektor( dNdX, dNdX);


		double** Y = new double*[4];
		for (int i = 0; i < 4; i++) {
			Y[i] = new double[4];
		}

		Y = wektorXwektor(dNdY, dNdY);


		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzH1[i][j] += this->k*(X[i][j] + Y[i][j])* this->Det_Jakobian[p];
			}
		}

	}//konczy sie petla po p calkowania w objetosci 

		//2 czesc macierzy H
		
		//identyczna tab f kszta³tu jak ta z local element
		double** N = new double*[8];
		for (int i = 0; i < 8; i++) {
			N[i] = new double[4];
		}

		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] = 0;
			}
		}


		//sprawdzanie warunków brzegowych
		if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[0] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[1] -1].status == 1) {	//warunek brzegowy spe³niony 
			N[0][0] = this->lE->N_P[0][0];
			N[0][1] = this->lE->N_P[0][1];
			N[0][2] = this->lE->N_P[0][2];
			N[0][3] = this->lE->N_P[0][3];

			N[1][0] = this->lE->N_P[1][0];
			N[1][1] = this->lE->N_P[1][1];
			N[1][2] = this->lE->N_P[1][2];
			N[1][3] = this->lE->N_P[1][3];
		}

		if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[1] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[2] -1].status == 1) {	//warunek brzegowy spe³niony 
			N[2][0] = this->lE->N_P[2][0];
			N[2][1] = this->lE->N_P[2][1];
			N[2][2] = this->lE->N_P[2][2];
			N[2][3] = this->lE->N_P[2][3];

			N[3][0] = this->lE->N_P[3][0];
			N[3][1] = this->lE->N_P[3][1];
			N[3][2] = this->lE->N_P[3][2];
			N[3][3] = this->lE->N_P[3][3];
		}

		if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[2] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[3] -1].status == 1) {	//warunek brzegowy spe³niony 
			N[4][0] = this->lE->N_P[4][0];
			N[4][1] = this->lE->N_P[4][1];
			N[4][2] = this->lE->N_P[4][2];
			N[4][3] = this->lE->N_P[4][3];

			N[5][0] = this->lE->N_P[5][0];
			N[5][1] = this->lE->N_P[5][1];
			N[5][2] = this->lE->N_P[5][2];
			N[5][3] = this->lE->N_P[5][3];
		}

		if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[3] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[0] -1].status == 1) {	//warunek brzegowy spe³niony 
			N[6][0] = this->lE->N_P[6][0];
			N[6][1] = this->lE->N_P[6][1];
			N[6][2] = this->lE->N_P[6][2];
			N[6][3] = this->lE->N_P[6][3];

			N[7][0] = this->lE->N_P[7][0];
			N[7][1] = this->lE->N_P[7][1];
			N[7][2] = this->lE->N_P[7][2];
			N[7][3] = this->lE->N_P[7][3];
		}

		double** NT = new double*[4];
		for (int j = 0; j < 4; j++) {
			NT[j] = new double[4];
		}
		for (int c = 0; c < 4; c++) {
			for (int d = 0; d < 4; d++) {
				NT[c][d] = 0;
			}
		}
		
		
		double detJ = 0; 

		for (int p = 0; p < 8; p++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					NT = wektorXwektor(N[p], N[p]);
					if (p == 0 || p == 1 || p == 6 || p == 7) {
						detJ = (this->B / (this->nB -1)) / 2;
					}
					else {
						detJ = (this->A / (this->nA -1)) / 2;
					}
					
					macierzH2[i][j] += NT[i][j] * this->alfa*detJ;
				}
			}

			for (int g = 0; g < 4; g++) {
				for (int f = 0; f < 4; f++) {
					NT[g][f] = 0;
				}
			}
		}


		//3 czêœæ macierzy H
		
		macierzH3 = Oblicz_C();

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzH3[i][j] /= this->krok_czasowy;
			}
		}
	

	//oblicz jakobian 
	//dosta³am 4x2  tab[0][0,1,2,3] po x 
	//tak samo po y
	//sumuje tab po x i y 
	//razy k i razy wyzn_jakobianu[p]
	//tablice = dodaje to co wysz³o} 


		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] = (macierzH1[i][j] + macierzH2[i][j] + macierzH3[i][j]);
			}
		}
		
		return H;
}

double* Grid::Oblicz_LP(int id_elem) {

	//macierz P
	double* P = new double[4];
	for (int i = 0; i < 4; i++) {
		P[i] = 0;
	}

	//identyczna tab f kszta³tu jak ta z local element
	double** N = new double*[8];
	for (int i = 0; i < 8; i++) {
		N[i] = new double[4];
	}

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 4; j++) {
			N[i][j] = 0;
		}
	}

	//sprawdzanie warunków brzegowych
	if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[0]-1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[1]-1].status == 1) {	//warunek brzegowy spe³niony 
		N[0][0] = this->lE->N_P[0][0];
		N[0][1] = this->lE->N_P[0][1];
		N[0][2] = this->lE->N_P[0][2];
		N[0][3] = this->lE->N_P[0][3];

		N[1][0] = this->lE->N_P[1][0];
		N[1][1] = this->lE->N_P[1][1];
		N[1][2] = this->lE->N_P[1][2];
		N[1][3] = this->lE->N_P[1][3];
	}

	if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[1] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[2]-1].status == 1) {	//warunek brzegowy spe³niony 
		N[2][0] = this->lE->N_P[2][0];
		N[2][1] = this->lE->N_P[2][1];
		N[2][2] = this->lE->N_P[2][2];
		N[2][3] = this->lE->N_P[2][3];

		N[3][0] = this->lE->N_P[3][0];
		N[3][1] = this->lE->N_P[3][1];
		N[3][2] = this->lE->N_P[3][2];
		N[3][3] = this->lE->N_P[3][3];
	}

	if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[2] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[3] -1].status == 1) {	//warunek brzegowy spe³niony 
		N[4][0] = this->lE->N_P[4][0];
		N[4][1] = this->lE->N_P[4][1];
		N[4][2] = this->lE->N_P[4][2];
		N[4][3] = this->lE->N_P[4][3];

		N[5][0] = this->lE->N_P[5][0];
		N[5][1] = this->lE->N_P[5][1];
		N[5][2] = this->lE->N_P[5][2];
		N[5][3] = this->lE->N_P[5][3];
	}

	if (this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[3] -1].status == 1 && this->tab_nodes[this->tab_elem[id_elem -1].id_nodes[0] -1].status == 1) {	//warunek brzegowy spe³niony 
		N[6][0] = this->lE->N_P[6][0];
		N[6][1] = this->lE->N_P[6][1];
		N[6][2] = this->lE->N_P[6][2];
		N[6][3] = this->lE->N_P[6][3];

		N[7][0] = this->lE->N_P[7][0];
		N[7][1] = this->lE->N_P[7][1];
		N[7][2] = this->lE->N_P[7][2];
		N[7][3] = this->lE->N_P[7][3];
	}

	double* P1 = new double[4];
	for (int i = 0; i < 4; i++) {
		P1[i] = 0;
	}

	double detJ = 0;
	for (int p = 0; p < 8; p++) {
		for (int n = 0; n < 4; n++) {
			if (p == 0 || p == 1 || p == 6 || p == 7) {
				detJ = (this->B / (this->nB - 1)) / 2;
			}
			else {
				detJ = (this->A / (this->nA - 1)) / 2;
			}
			P1[n] += (N[p][n] * this->alfa*this->t_otoczenia*detJ);
		}
	}


	//2 czêœæ macierzy P
	double** C = new double*[4];
	for (int i = 0; i < 4; i++) {
		C[i] = new double[4];
	}

	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			C[a][b] = 0;
		}
	}

	double** tempC = new double*[4];
	for (int i = 0; i < 4; i++) {
		tempC[i] = new double[4];
	}

	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			tempC[a][b] = 0;
		}
	}


	double* P2 = new double[4];
	double* T = new double[4];	//temp w i-tym p ca³kowania
	for (int i = 0; i < 4; i++) {
		T[i] = 0;
		P2[i] = 0;
	}


	for (int p = 0; p < 4; p++) {
		T[p] = (this->tab_nodes[this->tab_elem[id_elem - 1].id_nodes[0] - 1].temperature * this->lE->N_V[p][0]) + (this->tab_nodes[this->tab_elem[id_elem - 1].id_nodes[1] - 1].temperature * this->lE->N_V[p][1]) + (this->tab_nodes[this->tab_elem[id_elem - 1].id_nodes[2] - 1].temperature * this->lE->N_V[p][2]) + (this->tab_nodes[this->tab_elem[id_elem - 1].id_nodes[3] - 1].temperature * this->lE->N_V[p][3]);

	}


	for (int p = 0; p < 4; p++) {		//petla po p ca³kowania

		double* N = new double[4];
		for (int i = 0; i < 4; i++) {
			N[i] = 0;
		}

		N[0] = lE->N_V[p][0];
		N[1] = lE->N_V[p][1];
		N[2] = lE->N_V[p][2];
		N[3] = lE->N_V[p][3];

		tempC = wektorXwektor(N, N);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				C[i][j] = (tempC[i][j] * this->ro*this->cp*this->Det_Jakobian[p])/ this->krok_czasowy;
			}
		}


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			P2[i] += (T[p] * C[i][j]);
		}
	}


		for (int a = 0; a < 4; a++) {
			for (int b = 0; b < 4; b++) {
				tempC[a][b] = 0;
				C[a][b] = 0;
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		P[i] += (P1[i] + P2[i]);
	}

	return P;
}



double** Grid::Oblicz_C() {

	double** C = new double*[4];
	for (int i = 0; i < 4; i++) {
		C[i] = new double[4];
	}

	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			C[a][b] = 0;
		}
	}

	double** tempC = new double*[4];
	for (int i = 0; i < 4; i++) {
		tempC[i] = new double[4];
	}

	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			tempC[a][b] = 0;
		}
	}

	for (int p = 0; p < 4; p++) {		//petla po p ca³kowania
	
		double* N = new double[4];
		for (int i = 0; i < 4; i++) {
			N[i] = 0;
		}
	
		N[0] = lE->N_V[p][0];
		N[1] = lE->N_V[p][1];
		N[2] = lE->N_V[p][2];
		N[3] = lE->N_V[p][3];

		tempC = wektorXwektor(N, N);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				C[i][j] += tempC[i][j]*this->ro*this->cp*this->Det_Jakobian[p];
			}
		}
		
			for (int a = 0; a < 4; a++) {
				for (int b = 0; b < 4; b++) {
					tempC[a][b] = 0;
				}
			}
	}

	return C;
}

void Grid::Utworz_wezly(){//Node* tab_nodes) {

	//tab_nodes = new Node[this->num_nodes];
	//for (int i = 0; i < this->num_nodes; i++) {
	int i = 0;
		for (int n = 0; n < nB; n++) {
			for (int l = 0; l < nA; l++) {
				tab_nodes[i].id = i+1;
				tab_nodes[i].x = n * (B / (nB - 1));
				tab_nodes[i].y = l * (A / (nA - 1));
				tab_nodes[i].temperature = this->t0;
				if (n != 0 && n != (nB - 1) && l != 0 && l != (nA - 1)) {
					tab_nodes[i].status = false;
				}
				else {
					tab_nodes[i].status = true;
				}
				if (i < 16) {
					i++;
				}
				
			}
		}		
	//}		
		//for (int i = 0; i < num_nodes; i++)
		//	cout << tab_nodes[i].status << endl;

}

void Grid::Utworz_elementy(){//Element* tab_elem) {
	
	//tab_elem = new Element[this->num_elem];
	for (int i = 0; i < this->num_elem; i++) {
		
		tab_elem[i].id_element = i+1 ;//i+1

		if (i == 0 || i == 1 || i == 2) {
			tab_elem[i].id_nodes[0] = i + 1;
			tab_elem[i].id_nodes[1] = i + 5;
			tab_elem[i].id_nodes[2] = i + 6;
			tab_elem[i].id_nodes[3] = i + 2;
		}
		else if (i == 3 || i == 4 || i == 5) {
			tab_elem[i].id_nodes[0] = i + 2;
			tab_elem[i].id_nodes[1] = i + 6;
			tab_elem[i].id_nodes[2] = i + 7;
			tab_elem[i].id_nodes[3] = i + 3;
		}
		else {
			tab_elem[i].id_nodes[0] = i + 3;
			tab_elem[i].id_nodes[1] = i + 7;
			tab_elem[i].id_nodes[2] = i + 8;
			tab_elem[i].id_nodes[3] = i + 4;
		}
		for (int j = 0; j < 4; j++) {
			tab_elem[i].id_walls[j] = j + 1;
		}
	}
	
}


double** Grid::Oblicz_Jakobian(int id_elementu, int p) {	//p konkretny punkt ca³kowania
	

	double** Jakobian = new double*[2];
	for (int j = 0; j < 2; j++) {
		Jakobian[j] = new double[2];
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			Jakobian[i][j] = 0;
		}
	}
	
	//int index=this->tab_elem[id_elementu].id_nodes[p];

		Jakobian[0][0] = ((this->lE->dNdKsi[p][0])*this->tab_nodes[this->tab_elem[id_elementu - 1].id_nodes[0] -1].x + (this->lE->dNdKsi[p][1])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[1] -1].x + (this->lE->dNdKsi[p][2])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[2] -1].x + (this->lE->dNdKsi[p][3])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[3] -1].x);
		Jakobian[0][1] = ((this->lE->dNdKsi[p][0])*this->tab_nodes[this->tab_elem[id_elementu - 1].id_nodes[0] -1].y + (this->lE->dNdKsi[p][1])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[1] -1].y + (this->lE->dNdKsi[p][2])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[2] -1].y + (this->lE->dNdKsi[p][3])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[3] -1].y);
		Jakobian[1][0] = ((this->lE->dNdEta[p][0])*this->tab_nodes[this->tab_elem[id_elementu - 1].id_nodes[0] -1].x + (this->lE->dNdEta[p][1])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[1] -1].x + (this->lE->dNdEta[p][2])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[2] -1].x + (this->lE->dNdEta[p][3])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[3] -1].x);
		Jakobian[1][1] = ((this->lE->dNdEta[p][0])*this->tab_nodes[this->tab_elem[id_elementu - 1].id_nodes[0] -1].y + (this->lE->dNdEta[p][1])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[1] -1].y + (this->lE->dNdEta[p][2])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[2] -1].y + (this->lE->dNdEta[p][3])*this->tab_nodes[this->tab_elem[id_elementu -1].id_nodes[3] -1].y);
	//pochodne wrzuciæ do Jakobianu, s¹ obliczone w local element	
	//wyznacznik
		Det_Jakobian[p] = (Jakobian[0][0] * Jakobian[1][1]) - (Jakobian[0][1] * Jakobian[1][0]);
	//wyzn_jakobian[p] = wyznacznik dla tego punktu

		double** dNdXY = new double*[2];
		for (int k = 0; k < 2; k++) {
			dNdXY[k] = new double[4];
		}

		double** Jakobian_Odwrotny = new double*[2];
		for (int j = 0; j < 2; j++) {
			Jakobian_Odwrotny[j] = new double[2];
		}
		

		Jakobian_Odwrotny[0][0] = Jakobian[1][1];
		Jakobian_Odwrotny[0][1] = -Jakobian[0][1];
		Jakobian_Odwrotny[1][0] = -Jakobian[1][0];
		Jakobian_Odwrotny[1][1] = Jakobian[0][0];

		Jakobian_Odwrotny[0][0] = (1 / Det_Jakobian[p])*Jakobian_Odwrotny[0][0];
		Jakobian_Odwrotny[0][1] = (1 / Det_Jakobian[p])*Jakobian_Odwrotny[0][1];
		Jakobian_Odwrotny[1][0] = (1 / Det_Jakobian[p])*Jakobian_Odwrotny[1][0];
		Jakobian_Odwrotny[1][1] = (1 / Det_Jakobian[p])*Jakobian_Odwrotny[1][1];

		double** dNdKsiEta = new double*[2];
		for (int k = 0; k < 2; k++) {
			dNdKsiEta[k] = new double[4];
		}
		dNdKsiEta[0][0] = this->lE->dNdKsi[p][0];
		dNdKsiEta[0][1] = this->lE->dNdKsi[p][1];
		dNdKsiEta[0][2] = this->lE->dNdKsi[p][2];
		dNdKsiEta[0][3] = this->lE->dNdKsi[p][3];
		dNdKsiEta[1][0] = this->lE->dNdEta[p][0];
		dNdKsiEta[1][1] = this->lE->dNdEta[p][1];
		dNdKsiEta[1][2] = this->lE->dNdEta[p][2];
		dNdKsiEta[1][3] = this->lE->dNdEta[p][3];


		dNdXY = oblicz_macierze(Jakobian_Odwrotny, dNdKsiEta);


	//zapisaæ obliczone wartoœci w tablicy dndx i dndy //2x4

			//zwraca tablice dndx dndy 2x4

		return dNdXY;
			
}


double** Grid::oblicz_macierze(double** A, double** B) {
	double** C = new double*[2];
	for (int j = 0; j < 2; j++) {
		C[j] = new double[4];
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = 0;
		}
	}

	for (int k = 0; k < 2; k++) {
		for (int l = 0; l < 4; l++) {
			for (int i = 0; i < 2; i++) {
				C[k][l] += (A[k][i] * B[i][l]);		
			}
		}
	}
	return C;
}

double** Grid::wektorXwektor(double* A, double* B) {
	double** C = new double*[4];
	for (int j = 0; j < 4; j++) {
		C[j] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = 0;
		}
	}
	
		for (int l = 0; l < 4; l++) {
			for (int i = 0; i < 4; i++) {
				C[l][i] += (A[l] * B[i]);		
			}
		}
	
	return C;
}

double** Grid::MxM(double** A, double** B) {
	double** C = new double*[4];
	for (int j = 0; j < 4; j++) {
		C[j] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = 0;
		}
	}

	for (int k = 0; k < 4; k++) {
		for (int l = 0; l < 4; l++) {
			for (int i = 0; i < 8; i++) {
				C[k][l] += (A[k][i] * B[i][l]);
			}
		}
	}
	return C;
}