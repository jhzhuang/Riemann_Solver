#include <fstream>
#include "riemann_solver.h"

using namespace std;

void Solve(Discrete_Records &result, Status &s1, Status &s2);
void PrintResult(Discrete_Records &result);

int main(int argc, char **argv) {

	Status s1, s2;
	Discrete_Records result;
	

	s1.pressure = 1.0;
	s1.density = 1.0;
	s1.velocity = 0.0;

	s2.pressure = 0.1;
	s2.density = 0.125;
	s2.velocity = 0.0;

	result.records = new Status[201];
	result.time = 0.14;
	result.size = 201;
	result.dx = 2.0 / (201 - 1);

	Solve(result, s1, s2);

	PrintResult(result);

	delete[] result.records;

	return 0;
}

void Solve(Discrete_Records &result,Status &s1, Status &s2) {
	Riemann_Solver solver(s1,s2);
	solver.Solver();
	solver.Get_Result(result);
}

void PrintResult(Discrete_Records &result) {

	ofstream file("output.dat");

	file << "VARIABLES = \"X\", \"Pressure\", \"Density\", \"Velocity\"" << '\n';
	file << "ZONE I=201, DATAPACKING=BLOCK" << '\n';

	file.precision(10);
	file.setf(ios::fixed, ios::floatfield);

	for (int i = 0; i < result.size; i++)
	{
		file << i * result.dx - 1.0 << '\t';
	}
	file << '\n';
	for (int i = 0; i < result.size; i++)
	{
		file << result.records[i].pressure << '\t';
	}
	file << '\n';
	for (int i = 0; i < result.size; i++)
	{
		file << result.records[i].density << '\t';
	}
	file << '\n';
	for (int i = 0; i < result.size; i++)
	{
		file << result.records[i].velocity << '\t';
	}
	file << endl;
}
