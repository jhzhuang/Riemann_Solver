#include "pressure_solver.h"
#include "pressure_function.h"
#include "newton_solver.h"

double Modified_PF::Eval(double p) {
	return Pressure_Function(p, status_1, status_2) - (status_1.velocity - status_2.velocity);
}

double Pressure_Solver::Solver(double guess) {
	Newton_Solver<Modified_PF> p_solver(functor);
	return p_solver.Solver(guess);
}