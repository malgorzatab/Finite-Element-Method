#include "Local_Element.h"
#include <math.h>

#define ptk1_gauss 1/sqrt(3);
#define ptk2_gauss -1/sqrt(3);

Local_Element::Local_Element() {
	for (int i = 0; i < 4; i++) {
		this->ksi[i] = 0;
		this->eta[i] = 0;
	}

	for (int j = 0; j < 8; j++) {
		this->ksi_pow[j] = 0;
		this->eta_pow[j] = 0;
	}

	for (int k = 0; k < 4; k++) {
		for (int l = 0; l < 4; l++) {
			this->N_V[k][l] = 0;
			this->dNdKsi[k][l] = 0;
			this->dNdEta[k][l] = 0;

		}
	}

	for (int k = 0; k < 8; k++) {
		for (int l = 0; l < 4; l++) {
			this->N_P[k][l] = 0;
		}
	}
	
	

	setKsi();
	setEta();
	obliczN_V();
	obliczN_P();
	obliczdNdKsi();
	obliczdNdta();

}
//dokoñczyæ 
Local_Element::Local_Element(double* ksi, double*eta, double* ksi_pow, double* eta_pow, double**N_V, double** N_P, double** dNdKsi, double** dNdEta)
{
	for (int i = 0; i < 4; i++) {
		this->ksi[i] = ksi[i];
		eta[i] = 0;
	}
	
	for (int j = 0; j < 8; j++) {
		ksi_pow[j] = 0;
		eta_pow[j] = 0;
	}

	for (int k = 0; k < 4; k++) {
		for (int l = 0; l < 4; l++) {
			N_V[k][l] = 0;
			N_P[k][l] = 0;
			dNdKsi[k][l] = 0;
			dNdEta[k][l] = 0;

		}
	}

	setKsi();
	setEta();
	obliczN_V( );
	obliczN_P();
	obliczdNdKsi();
	obliczdNdta();
}

Local_Element::~Local_Element(){}

void Local_Element::setKsi() {
	this->ksi[0] = ptk2_gauss;
	this->ksi[1] = ptk1_gauss;
	this->ksi[2] = ptk1_gauss;
	this->ksi[3] = ptk2_gauss;
}
void Local_Element::setEta() {
	this->eta[0] = ptk2_gauss;
	this->eta[1] = ptk2_gauss;
	this->eta[2] = ptk1_gauss;
	this->eta[3] = ptk1_gauss;
}

void Local_Element::obliczN_V() {	
	//setKsi();
	//setEta();

	for (int i = 0; i < 4; i++) {		
		this->N_V[i][0] = 0.25*(1 - this->ksi[i])*(1 - this->eta[i]);
		this->N_V[i][1] = 0.25*(1 + this->ksi[i])*(1 - this->eta[i]);
		this->N_V[i][2] = 0.25*(1 + this->ksi[i])*(1 + this->eta[i]);
		this->N_V[i][3] = 0.25*(1 - this->ksi[i])*(1 + this->eta[i]);
	}		
}

void Local_Element::obliczN_P() {
	this->ksi_pow[0] = ptk2_gauss;
	this->ksi_pow[1] = ptk1_gauss;
	this->ksi_pow[2] = 1;
	this->ksi_pow[3] = 1;
	this->ksi_pow[4] = ptk1_gauss;
	this->ksi_pow[5] = ptk2_gauss;
	this->ksi_pow[6] = -1;
	this->ksi_pow[7] = -1;

	this->eta_pow[0] = -1;
	this->eta_pow[1] = -1;
	this->eta_pow[2] = ptk2_gauss;
	this->eta_pow[3] = ptk1_gauss;
	this->eta_pow[4] = 1;
	this->eta_pow[5] = 1;
	this->eta_pow[6] = ptk1_gauss;
	this->eta_pow[7] = ptk2_gauss;

	for (int i = 0; i < 8; i++) {
		this->N_P[i][0] = 0.25*(1 - this->ksi_pow[i])*(1 - this->eta_pow[i]);
		this->N_P[i][1] = 0.25*(1 + this->ksi_pow[i])*(1 - this->eta_pow[i]);
		this->N_P[i][2] = 0.25*(1 + this->ksi_pow[i])*(1 + this->eta_pow[i]);
		this->N_P[i][3] = 0.25*(1 - this->ksi_pow[i])*(1 + this->eta_pow[i]);

	}
}

void Local_Element::obliczdNdKsi() {
	//setEta();
	for (int i = 0; i < 4; i++) {
		this->dNdKsi[i][0] = -0.25*(1 - this->eta[i]);
		this->dNdKsi[i][1] = 0.25*(1 - this->eta[i]);
		this->dNdKsi[i][2] = 0.25*(1 + this->eta[i]);
		this->dNdKsi[i][3] = -0.25*(1 + this->eta[i]);
	}

}
void Local_Element::obliczdNdta() {
	//setKsi();
	for (int i = 0; i < 4; i++) {
		this->dNdEta[i][0] = -0.25*(1 - this->ksi[i]);
		this->dNdEta[i][1] = -0.25*(1 + this->ksi[i]);
		this->dNdEta[i][2] = 0.25*(1 + this->ksi[i]);
		this->dNdEta[i][3] = 0.25*(1 - this->ksi[i]);
	}
}