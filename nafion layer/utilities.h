#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED
#include "mesh.h"
#include "potential_signal.h"
#include "thermodynamics.h"
#include <memory>

void initialize(mesh& membrane, mesh& solution, nernst_equation& therm, const double& Eq, const double& Cm_bulk);
void membrane_matrixA(mesh& membrane, const potential_signal& Esignal, const double D_Fcm, const double D_Fcpm, const double K_dis);
void solution_matrixA(mesh& solution, const potential_signal& Esignal, const mesh& membrane, const double D_Fcw, const double D_Fcpw, const double K_dis);
void membrane_arrayb(mesh& membrane);
void solution_arrayb(mesh& solution);
double solver(const int ti, const double& Eq, mesh& membrane, mesh& solution, nernst_equation& therm, const potential_signal& Esignal);
#endif