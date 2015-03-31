#ifndef _RIEMANNSOLVER_PRESSURE_FUNCTION_H_
#define _RIEMANNSOLVER_PRESSURE_FUNCTION_H_

#include "gas.h"
#include <math.h>

#define Gamma 1.4

/* Unilateral Pressure Function */
inline double Unilateral_PF(double p_star, double p, double density) {
	if (p_star > p)
	{
		return (p_star - p) / (density * (sqrt((p / density) * (0.5 * (Gamma + 1) * (p_star / p) + 0.5 * (Gamma - 1)))));
	}
	else
	{
		return (2.0 / (Gamma - 1)) * sqrt(Gamma * p / density) * (pow((p_star / p), (Gamma - 1) / (2 * Gamma)) - 1);
	}
}

#ifdef _RIEMANNSOLVER_PRESSURE_SOLVER_H_

double Pressure_Function(double p_star, Status &status_1, Status &status_2) {
	return Unilateral_PF(p_star, status_1.pressure, status_1.density) + Unilateral_PF(p_star, status_2.pressure, status_2.density);
}

#endif

#endif