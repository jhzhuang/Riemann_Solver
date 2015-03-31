#ifndef _RIEMANNSOLVER_NEWTON_SOLVER_H_
#define _RIEMANNSOLVER_NEWTON_SOLVER_H_

#include <math.h>

#define Infinitestimal 1.0e-7

#define Confine_Infinitestimal(X) (Infinitestimal > abs(X) * Infinitestimal ? Infinitestimal : abs(X) * Infinitestimal )

template<typename Functor>
class Newton_Solver {
public:
	Newton_Solver(Functor &fun) : functor(fun) {}
	~Newton_Solver() {}
public:
	double Solver(double guess);
private:
	Functor &functor;
};

template<typename Functor>
double Newton_Solver<Functor>::Solver(double guess) {

	double dx = Confine_Infinitestimal(guess);
	double x = guess;

	while (abs(functor.Eval(x)) > Infinitestimal)
	{
		x = x - functor.Eval(x) * (2.0 * dx / (functor.Eval(x + dx) - functor.Eval(x - dx)));
		dx = Confine_Infinitestimal(x);
	}

	return x;
}

#endif