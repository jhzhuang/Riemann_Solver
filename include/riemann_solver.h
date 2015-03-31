#ifndef _RIEMANNSOLVER_RIEMANN_SOLVER_H_
#define _RIEMANNSOLVER_RIEMANN_SOLVER_H_

#include "physical_objects.h"
#include "field.h"

typedef struct _Discrete_Records {
	Status *records;
	double time;
	double dx;
	int size;    // varialble size should be odd and no less than 3
} Discrete_Records;

class Riemann_Solver {
public:
	Riemann_Solver(Status &status_1, Status &status_2);
	~Riemann_Solver() {}
public:
	void Solver();
	void Get_Result(Discrete_Records &result);    //using function Get_Result before call function Solver
private:
	Status_Container container;
	Field field;
};

#endif