#include <stdio.h>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Element.h"
using namespace std;

Element::Element() {
	this->id_element = 0;
	for (int i = 0; i < 4; i++) {
		this->id_nodes[i] = 0;
		this->id_walls[i] = 0;
		this->lP[i] = 0;
	}

	this->lH = new double*[4];
	for (int i = 0; i < 4; i++) {
		this->lH[i] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->lH[i][j] = 0;
		}
	}
}

Element::Element(int id, int* walls, int* id_nodes, double** lH, double* lP) {
	this->id_element = id;
	for (int i = 0; i < 4; i++) {
		this->id_nodes[i] = id_nodes[i];
		this->id_walls[i] = walls[i];
		this->lP[i] = 0;
	}

	this->lH = new double*[4];
	for (int i = 0; i < 4; i++) {
		this->lH[i] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->lH[i][j] = 0;		//lH[i][j]
		}
	}
}

Element::~Element() {}


void Element::setlH(double** lH) {

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->lH[i][j] = lH[i][j];
		}
	}
}

void Element::setlP(double* lP) {
		for (int i = 0; i < 4; i++) {
			this->lP[i] = lP[i];
		}
	}


//robi to samo co setlH xd
void Element::set_localH(double** lH) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->lH[i][j] = lH[i][j];
		}
	}
}

void Element::set_localP(double* lP) {
	for (int i = 0; i < 4; i++) {
		this->lP[i] = lP[i];
	}
}

void Element::setId_element(int id_element) {
	this->id_element = id_element;
}

void Element::setId_walls(int* id_walls) {
	for (int i = 0; i < 4; i++) {
		this->id_walls[i] = i + 1;
	}
}

void Element::setId_nodes(int* id_nodes) {
	for (int i = 0; i < 4; i++) {
		this->id_nodes[i] = i + 1;
	}
}