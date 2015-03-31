#ifndef _RIEMANNSOLVER_PRESSURE_SOLVER_H_
#define _RIEMANNSOLVER_PRESSURE_SOLVER_H_

#include "gas.h"

/* Modified Pressure Function */
class Modified_PF {
public:
	Modified_PF(Status &s1, Status &s2) : status_1(s1), status_2(s2) {}
	~Modified_PF() {}
public:
	double Eval(double p);    //interface for Newton solver
private:
	Status &status_1;
	Status &status_2;
};

class Pressure_Solver {
public:
	Pressure_Solver(Status &s1, Status &s2) : functor(s1, s2) {}
	~Pressure_Solver() {}
public:
	double Solver(double guess);
private:
	Modified_PF functor;
};

#endif