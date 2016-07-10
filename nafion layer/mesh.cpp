#include "mesh.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

mesh::mesh(int fm, double fdr0, double fdr, int fn, double fdz0, double fdz) :
m(fm), dr0(fdr0), dr(fdr), n(fn), dz0(fdz0), dz(fdz)
{
	R = new double *[m];
	RR = new double *[m];
	Z = new double *[m];
	ZZ = new double *[m];
	Crdn = new double *[m];
	Crdo = new double *[m];
	Coxn = new double *[m];
	Coxo = new double *[m];
	arrayb = new double *[m+2];// add n elements at the beginning and the end of the list respectively
	arrayb2 = new double *[m+2];
	matrixA = new double **[m+2];
	matrixA2 = new double **[m+2];
	x = new double *[m + 2];
	xk = new double *[m + 2];
	xo = new double *[m + 2];
	xp = new double *[m + 2];
	xkp = new double *[m + 2];
	xpo = new double *[m + 2];


	for (int i = 0; i < m; ++i){
		R[i] = new double [n];
		RR[i] = new double [n];
		Z[i] = new double [n];
		ZZ[i] = new double [n];
		Crdn[i] = new double [n];
		Crdo[i] = new double [n];
		Coxn[i] = new double [n];
		Coxo[i] = new double [n];
	}

	for (int i = 0; i < m + 2; ++i){
		arrayb[i] = new double[n + 2];
		arrayb2[i] = new double[n + 2];
		matrixA[i] = new double *[n + 2];
		matrixA2[i] = new double *[n + 2];
		x[i] = new double[n + 2];
		xk[i] = new double[n + 2];
		xo[i] = new double[n + 2];
		xp[i] = new double[n + 2];
		xkp[i] = new double[n + 2];
		xpo[i] = new double[n + 2];

		for (int j = 0; j < n + 2; ++j){
			arrayb[i][j] = 0;
			arrayb2[i][j] = 0;
			x[i][j] = 0;
			xk[i][j] = 0;
			xo[i][j] = 0;
			xp[i][j] = 0;
			xkp[i][j] = 0;
			xpo[i][j] = 0;
		}
	}

	for (int i = 0; i < m + 2; ++i){
		for (int j = 0; j < n + 2; ++j){
			matrixA[i][j] = new double[5];
			matrixA2[i][j] = new double[5];
			for (int k = 0; k < 5; ++k){
				matrixA[i][j][k] = 0;
				matrixA2[i][j][k] = 0;
			}
		}
		
	}

	// generate mesh grid
	mesh_grid();

	// cacluate late geometric coefficient for diffusion
	alfa1 = new double [m];
	alfa2 = new double [m];
	alfa3 = new double [m];
	for (int i = 0; i < m; ++i){
		if (i == 0){
			alfa1[i] = 2 / (3 * dr0 + dr)*(1 / dr0 - 1 / RR[i][0]);
			alfa2[i] = 2 / (3 * dr0 + dr)*(2 / (dr0 + dr) + 1 / RR[i][0]);
			alfa3[i] = -2 / dr0 / (dr0 + dr);
		}
		else if (i == 1){
			alfa1[i] = 2 / (dr0 + 3 * dr)*(2 / (dr0 + dr) - 1 / RR[i][0]);
			alfa2[i] = 2 / (dr0 + 3 * dr)*(2 / (dr0 + dr) + 1 / RR[i][0]);
			alfa3[i] = -2 / (dr0 + dr) / dr;
		}
		else{
			alfa1[i] = 0.5 / dr*(1 / dr - 1 / RR[i][0]);
			alfa2[i] = 0.5 / dr*(1 / dr + 1 / RR[i][0]);
			alfa3[i] = -1 / dr / dr;
		}
	}
	beta1 = new double[n];
	beta2 = new double[n];
	beta3 = new double[n];
	for (int i = 0; i < n; ++i){
		if (i == 0){
			beta1[i] = 2 / dz0 / (3 * dz0 + dz);
			beta2[i] = 4 / (dz0 + dz) / (3 * dz0 + dz);
			beta3[i] = -2 / dz0 / (dz0 + dz);
		}
		else if (i == 1){
			beta1[i] = 4 / (dz0 + dz) / (dz0 + 3 * dz);
			beta2[i] = 2 / dz / (dz0 + 3 * dz);
			beta3[i] = -2 / dz0 / (dz0 + dz);
		}
		else{
			beta1[i] = 0.5 / dz / dz;
			beta2[i] = 0.5 / dz / dz;
			beta3[i] = -1 / dz / dz;
		}
	}
}

mesh::~mesh()
{
	for (int i = 0; i < m; ++i){
		delete R[i];
		delete Z[i];
		delete RR[i];
		delete ZZ[i];
		delete Crdn[i];
		delete Crdo[i];
		delete Coxn[i];
		delete Coxo[i];
	}

	for (int i = 0; i < m + 2; ++i){
		for (int j = 0; j < n; ++j){
			delete matrixA[i][j];
			delete matrixA2[i][j];
		}
	}

	for (int i = 0; i < m + 2; ++i){
		delete matrixA[i];
		delete matrixA2[i];
		delete arrayb[i];
		delete arrayb2[i];
		delete x[i];
		delete xk[i];
		delete xkp[i];
		delete xo[i];
		delete xp[i];
		delete xpo[i];
	}

	delete R;
	delete Z;
	delete RR;
	delete ZZ;
	delete Crdn;
	delete Crdo;
	delete Coxn;
	delete Coxo;
	delete matrixA;
	delete matrixA2;
	delete arrayb;
	delete arrayb2;
	delete alfa1;
	delete alfa2;
	delete alfa3;
	delete beta1;
	delete beta2;
	delete beta3;
	delete x;
	delete xk;
	delete xkp;
	delete xo;
	delete xp;
	delete xpo;
}

void mesh::mesh_grid()
{
	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j){
			if (i != 0 && j != 0){
				R[i][j] = (i - 1)*dr + dr0;
				Z[i][j] = (j - 1)*dz + dz0;
				RR[i][j] = R[i][j] + 0.5*dr;
				ZZ[i][j] = Z[i][j] + 0.5*dz;
				Crdn[i][j] = 0;
				Crdo[i][j] = 0;
				Coxn[i][j] = 0;
				Coxo[i][j] = 0;
			}
			else if (i == 0 && j != 0){
				R[i][j] = 0;
				Z[i][j] = (j - 1)*dz + dz0;
				RR[i][i] = 0.5*dr0;
				ZZ[i][j] = Z[i][j] + 0.5*dz;
				Crdn[i][j] = 0;
				Crdo[i][j] = 0;
				Coxn[i][j] = 0;
				Coxo[i][j] = 0;
			}
			else if (i != 0 && j == 0){
				R[i][j] = (i - 1)*dr + dr0;
				Z[i][j] = 0;
				RR[i][j] = R[i][j] + 0.5*dr;
				ZZ[i][j] = 0.5*dz0;
				Coxn[i][j] = 0;
				Coxo[i][j] = 0;
			}
			else if (i == 0 && j == 0){
				R[i][j] = 0;
				Z[i][j] = 0;
				RR[i][i] = 0.5*dr0;
				ZZ[i][j] = 0.5*dz0;
				Crdn[i][j] = 0;
				Crdo[i][j] = 0;
				Coxn[i][j] = 0;
				Coxo[i][j] = 0;
			}
		}
	}
}

void mesh::Ca2Cb(string Ca_name, string Cb_name)
{
	double **Ca, **Cb;
	if (Ca_name == "Crdn"){
		Ca = Crdn;
	}
	else if (Ca_name == "Crdo"){
		Ca = Crdo;
	}
	else if (Ca_name == "Coxn"){
		Ca = Coxn;
	}
	else if (Ca_name == "Coxo"){
		Ca = Coxo;
	}
	else{
		std::cout << "Ca2Cb cannot find Ca: " << Ca_name << endl;
		exit(EXIT_FAILURE);

	}

	if (Cb_name == "Crdn"){
		Cb = Crdn;
	}
	else if (Cb_name == "Crdo"){
		Cb = Crdo;
	}
	else if (Cb_name == "Coxn"){
		Cb = Coxn;
	}
	else if (Cb_name == "Coxo"){
		Cb = Coxo;
	}
	else{
		std::cout << "Ca2Cb cannot find Cb: " << Cb_name << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j){
			Cb[i][j] = Ca[i][j];
		}
	}
}

void mesh::print_mesh(string ex_file_name)
{
	ofstream myfile1(ex_file_name + " mesh-R.txt");
	ofstream myfile2(ex_file_name + " mesh-Z.txt");
	ofstream myfile3(ex_file_name + " mesh parameters.txt");

	if (myfile1.is_open() && myfile2.is_open()){
		for (int j = 0; j < n; ++j){
			for (int i = 0; i < m; ++i){
				myfile1 << R[i][j] << " ";
				myfile2 << Z[i][j] << " ";
			}
			myfile1 << endl;
			myfile2 << endl;
		}
		myfile1.close();
		myfile2.close();
	}
	else{
		cout << " unable to open mesh data file" << endl;
		exit(EXIT_FAILURE);
	}

	myfile3 << "number of R nodes: " << m
		<< "\ndr0, cm: " << dr0
		<< "\ndr, cm: " << dr
		<< "number of z nodes: " << n
		<< "\ndz0, cm: " << dz0
		<< "\ndz, cm: " << dz;
}

void mesh::print_concentration(string target_concentration_name, string ex_file_name)
{
	double **C;
	if (target_concentration_name == "Crdn"){
		C = Crdn;
	}
	else if (target_concentration_name == "Crdo"){
		C = Crdo;
	}
	else if (target_concentration_name == "Coxn"){
		C = Coxn;
	}
	else if (target_concentration_name == "Coxo"){
		C = Coxo;
	}
	else{
		std::cout << "print_concentration cannot find target_concentration_name: " << target_concentration_name << endl;
		exit(EXIT_FAILURE);
	}

	ofstream myfile(ex_file_name + ".txt");
	if (myfile.is_open()){
		for (int j = 0; j < n; ++j){
			for (int i = 0; i < m; ++i){
				myfile << C[i][j] << " ";
			}
			myfile << endl;
		}
		myfile.close();
	}
	else{
		std::cout << "print_concentration cannot open " << target_concentration_name << ".txt" << endl;
	}
}