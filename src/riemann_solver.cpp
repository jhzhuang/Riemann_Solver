#include "riemann_solver.h"
#include "pressure_function.h"
#include "pressure_solver.h"

#define IS_SHOCK_WAVE(p_front,p_back) (p_front < p_back)
#define IS_VACCUM(p) (p <= 0.0)
#define IS_FILLED(p) (p > 0.0)
#define IS_LEFT(wave) ((wave.back - wave.front) == 1)

void Density_Solver(Status_Container &container, Physical_Object &wave);
void Discrete_Solver(Status_Container &container, Physical_Object &wave, Discrete_Records &record);

Riemann_Solver::Riemann_Solver(Status &status_1, Status &status_2) {

	container.push(status_1);
	container.push(status_1);
	container.push(status_2);
	container.push(status_2);

	field.left.front = 0;
	field.left.back = 1;
	field.right.front = 3;
	field.right.back = 2;
}

void Riemann_Solver::Solver() {

	double p_star;
	double u_star;
	Pressure_Solver p_solver(container.at(0),container.at(3));

	p_star = p_solver.Solver(0.5 * (container.at(0).pressure + container.at(3).pressure));
	if (IS_VACCUM(p_star))
	{
		p_star = 0.0;
	}

	u_star = 0.5 * (container.at(0).velocity + container.at(3).velocity + Unilateral_PF(p_star, container.at(3).pressure, container.at(3).density) - Unilateral_PF(p_star, container.at(0).pressure, container.at(0).density));
	
	container.at(1).pressure = p_star;
	container.at(1).velocity = u_star;
	Density_Solver(container, field.left);

	container.at(2).pressure = p_star;
	container.at(2).velocity = u_star;
	Density_Solver(container, field.right);
}

void Riemann_Solver::Get_Result(Discrete_Records &result) {
	Discrete_Solver(container, field.left, result);
	Discrete_Solver(container, field.right, result);
}

void Density_Solver(Status_Container &container, Physical_Object &wave) {
	if (IS_SHOCK_WAVE(container.at(wave.front).pressure, container.at(wave.back).pressure))
	{
		double p_ratio = container.at(wave.back).pressure / container.at(wave.front).pressure;
		double alpha = (Gamma + 1.0) / (Gamma - 1.0);
		container.at(wave.back).density = container.at(wave.front).density * ((alpha * p_ratio + 1.0) / (alpha + p_ratio));
		wave.type = shock_wave;
	}
	else if (IS_FILLED(container.at(wave.back).pressure))
	{
		double p_ratio = container.at(wave.back).pressure / container.at(wave.front).pressure;
		container.at(wave.back).density = container.at(wave.front).density * pow(p_ratio, 1.0 / Gamma);
		wave.type = expansion_wave;
	}
	else
	{
		container.at(wave.back).density = 0.0;
		wave.type = expansion_wave;
	}
}

void Discrete_Solver(Status_Container &container, Physical_Object &wave, Discrete_Records &record) {

	int i;
	int pesudo_i;

	double dx = record.dx;
	Status *discrete = record.records;
	int size = record.size;
	double time = record.time;

	int di = wave.back - wave.front;
	int left_end = - (size - 1) / 2;
	int right_end = (size - 1) / 2;

	if (IS_SHOCK_WAVE(container.at(wave.front).pressure, container.at(wave.back).pressure))
	{
		double u_shock = (container.at(wave.back).density * container.at(wave.back).velocity - container.at(wave.front).density * container.at(wave.front).velocity) / (container.at(wave.back).density - container.at(wave.front).density);
		int split_1 = (int)(u_shock * time / dx) * di;
		int split_2 = (int)(container.at(wave.back).velocity * time / dx) * di;

		for (pesudo_i = left_end, i = left_end * di + right_end; pesudo_i < split_1; pesudo_i++)
		{
			discrete[i] = container.at(wave.front);
			i = i + di;
		}

		for (pesudo_i = split_1; pesudo_i <= split_2; pesudo_i++)
		{
			discrete[i] = container.at(wave.back);
			i = i + di;
		}
	}
	else if (IS_FILLED(container.at(wave.back).pressure))
	{
		double u_wave_front = container.at(wave.front).velocity - di * sqrt(Gamma * container.at(wave.front).pressure / container.at(wave.front).density);
		double u_wave_back = container.at(wave.back).velocity - di * sqrt(Gamma * container.at(wave.back).pressure / container.at(wave.back).density);
		int split_1 = (int)(u_wave_front * time / dx) * di;
		int split_2 = (int)(u_wave_back * time / dx) * di;
		int split_3 = (int)(container.at(wave.back).velocity * time / dx) * di;

		for (pesudo_i = left_end, i = left_end * di + right_end; pesudo_i < split_1; pesudo_i++)
		{
			discrete[i] = container.at(wave.front);
			i = i + di;
		}

		for (pesudo_i = split_1; pesudo_i < split_2; pesudo_i++)
		{
			int xi = i + left_end;
			double c_front = sqrt(Gamma * container.at(wave.front).pressure / container.at(wave.front).density);
			double c = (Gamma - 1.0) / (Gamma + 1.0) * (container.at(wave.front).velocity - xi * dx / time) * di + (2.0 / (Gamma + 1.0)) * c_front;
			discrete[i].velocity = xi * dx / time + c * di;
			discrete[i].pressure = container.at(wave.front).pressure * pow(c / c_front, 2.0 * Gamma / (Gamma - 1.0));
			discrete[i].density = Gamma * discrete[i].pressure / (c * c);
			i = i + di;
		}

		for (pesudo_i = split_2; pesudo_i <= split_3; pesudo_i++)
		{
			discrete[i] = container.at(wave.back);
			i = i + di;
		}
	}
	else
	{
		double u_wave_front = container.at(wave.front).velocity - di * sqrt(Gamma * container.at(wave.front).pressure / container.at(wave.front).density);
		double u_wave_back = container.at(wave.front).velocity - di * 2.0 * sqrt(Gamma * container.at(wave.front).pressure / container.at(wave.front).density) / (Gamma + 1.0);
		int split_1 = (int)(u_wave_front * time / dx) * di;
		int split_2 = (int)(u_wave_back * time / dx) * di;
		int split_3 = (int)(container.at(wave.back).velocity * time / dx) * di;

		for (pesudo_i = left_end, i = left_end * di + right_end; pesudo_i < split_1; pesudo_i++)
		{
			discrete[i] = container.at(wave.front);
			i = i + di;
		}

		for (pesudo_i = split_1; pesudo_i < split_2; pesudo_i++)
		{
			int xi = i + left_end;
			double c_front = sqrt(Gamma * container.at(wave.front).pressure / container.at(wave.front).density);
			double c = (Gamma - 1.0) / (Gamma + 1.0) * (container.at(wave.front).velocity - xi * dx / time) * di + (2.0 / (Gamma + 1.0)) * c_front;
			discrete[i].velocity = xi * dx / time + c * di;
			discrete[i].pressure = container.at(wave.front).pressure * pow(c / c_front, 2.0 * Gamma / (Gamma - 1.0));
			discrete[i].density = Gamma * discrete[i].pressure / (c * c);
			i = i + di;
		}

		for (pesudo_i = split_2; pesudo_i <= split_3; pesudo_i++)
		{
			discrete[i].velocity = 0.0;
			discrete[i].pressure = 0.0;
			discrete[i].density = 0.0;
			i = i + di;
		}
	}
}