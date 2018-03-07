#include "GlobalData.h"
#include <stdio.h>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
using namespace std;

const double eps = 1e-12;

GlobalData::GlobalData() {
	this->num_nodes = 0;
	this->num_elem = 0;
	this->t0 = 0;
	this->max_czas = 0;
	this->time_step = 0;
	this->t_otoczenia = 0;
	this->alfa = 0;
	this->A = 0;
	this->B = 0;
	this->nA = 0;
	this->nB = 0;
	this->cp = 0;
	this->k = 0;
	this->ro = 0;

	this->global_H = new double*[16];
	for (int k = 0; k < 16; k++) {
		this->global_H[k] = new double[16];
	}
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			this->global_H[i][j] = 0;
		}
	}
	this->global_P = new double[16];
	for (int j = 0; j < 16; j++) {
		this->global_P[j] = 0;
	}

	this->temperatury = new double[16];
	for (int j = 0; j < 16; j++) {
		this->temperatury[j] = 0;
	}
	
}

//dokoñczyæ
GlobalData::GlobalData(int num_nodes, int num_elem, double** globalH, double* globalP, double* temp) {
	this->num_nodes = num_nodes;
	this->num_elem = num_elem;

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			this->global_H[i][j] = 0;
		}
	}

	for (int j = 0; j < 16; j++) {
		this->global_P[j] = 0;
	}

	for (int j = 0; j < 16; j++) {
		this->temperatury[j] = 0;
	}

}


void GlobalData::wczytaj() {
	int i = 0;
	string::size_type sz;
	string lines[14];
	fstream file;

	file.open("dane.txt", ios::in);

	if (file.good() == true) {
		while (!file.eof() && i < 14) {
			if (i < 14) {
				getline(file, lines[i]);
			}
			i++;
			file.flush();
		}

		file.close();
	}
	else {
		cout << "Nie mo¿na odczytaæ pliku!" << endl;
	}
	this->t0 = stod(lines[0], &sz);
	this->max_czas = stod(lines[1], &sz);
	this->time_step = stod(lines[2], &sz);
	this->t_otoczenia = stod(lines[3], &sz);
	this->alfa = stod(lines[4], &sz);
	this->A = stod(lines[5], &sz);
	this->B = stod(lines[6], &sz);
	this->nA = atoi(lines[7].c_str());
	this->nB = atoi(lines[8].c_str());
	this->cp = stod(lines[9], &sz);
	this->k = stod(lines[10], &sz);
	this->ro = stod(lines[11], &sz);
	this->num_elem = atoi(lines[12].c_str());
	this->num_nodes = atoi(lines[13].c_str());

}

void GlobalData::oblicz_GH(Grid grid){	
	grid.Ustaw_LH();
	for (int i = 0; i < 9; i++) {	
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				global_H[grid.tab_elem[i].id_nodes[j] -1][grid.tab_elem[i].id_nodes[k] -1 ] += grid.tab_elem[i].lH[j][k]; 
			}
		}
		
	}

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			cout << setw(3) << global_H[i][j] << "\t";
		}
		cout << endl;
	}
}

void GlobalData::oblicz_GP(Grid grid) {
	grid.Ustaw_LP();
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 4; j++) {
			global_P[grid.tab_elem[i].id_nodes[j] -1] += grid.tab_elem[i].lP[j];
		}
	}


	for (int j = 0; j < 16; j++) {
		cout << setw(3) << global_P[j] << "\t";
	}
}

void GlobalData::utworz_Siatke() {
	wczytaj();
	//Grid grid = Grid();
	Grid grid = Grid(this->A, this->B, this->nA, this->nB, this->num_nodes, this->num_elem, this->t0, this->t_otoczenia, this->max_czas, this->time_step, this->alfa, this->cp, this->k, this->ro);

	for (int i = 0; i < 16; i++) {
		this->temperatury[i] = this->t0;
	}

	int iter = 0;
	for (int t = 100; t <= max_czas; t += time_step) {
		cout << "Iteracja: " << iter << endl;
		cout << "Macierz H" << endl;
		cout << endl;
		oblicz_GH(grid);
		cout << endl;
		cout << "Macierz P" << endl;
		cout << endl;
		oblicz_GP(grid);
		cout << endl;
		this->temperatury = gauss(global_H, global_P);

		cout << "Temperatury: " << endl;
		for (int i = 0; i < 16; i++) {
			cout << "Time:		Temp:		" << endl;
			cout << t << "		" << this->temperatury[i] << endl;
		}

		for (int k = 0; k < 16; k++) {
			grid.tab_nodes[k].temperature = this->temperatury[k];
		}

		//wyzerowanie macierzy H i P
		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 16; j++) {
				this->global_H[i][j] = 0;
			}
		}
		
		for (int j = 0; j < 16; j++) {
			this->global_P[j] = 0;
			this->temperatury[j] = 0;
		}
iter++;
	}
	
}

double* GlobalData::gauss( double** H, double* P) {
	//int i, j, k;
	//double m, s;

	// eliminacja wspó³czynników

	/*for (i = 0; i < 16 - 1; i++)
	{
		for (j = i + 1; j < 16; j++)
		{
			if (fabs(H[i][i]) < eps) return false;
			m = -H[j][i] / H[i][i];
			for (k = i + 1; k <= 16; k++)
				H[j][k] += m * H[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = 16 - 1; i >= 0; i--)
	{
		s = H[i][16];
		for (j = 16 - 1; j >= i + 1; j--)
			s -= H[i][j] * P[j];
		if (fabs(H[i][i]) < eps) return false;
		P[i] = s / H[i][i];
	}*/

	double m, s;
	double eps = 1e-12;

	double** gG = new double *[num_nodes];
	for (int i = 0; i <= num_nodes; i++) {
		gG[i] = new double[num_nodes + 1];
	}

	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			gG[i][j] = H[i][j];
		}
		gG[i][num_nodes] = P[i];
	}

	for (int i = 0; i < num_nodes - 1; i++)
	{
		for (int j = i + 1; j < num_nodes; j++)
		{
			if (fabs(gG[i][i]) < eps) return false;
			m = -gG[j][i] / gG[i][i];
			for (int k = i + 1; k <= num_nodes; k++)
				gG[j][k] += m * gG[i][k];
		}
	}

	for (int i = num_nodes - 1; i >= 0; i--)
	{
		s = gG[i][num_nodes];
		for (int j = num_nodes - 1; j >= i + 1; j--)
			s -= gG[i][j] * this->temperatury[j];
		if (fabs(gG[i][i]) < eps) return false;
		this->temperatury[i] = s / gG[i][i];

	}
	return this->temperatury;
}