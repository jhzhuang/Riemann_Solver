#include <fstream>
#include "riemann_solver.h"

using namespace std;

void Solve(Discrete_Records &result, Status &s1, Status &s2);

int main(int argc, int argv) {

	ofstream file("output.dat");
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

	file << "VARIABLES = \"X\", \"Pressure\", \"Density\", \"Velocity\"" << endl;
	file << "ZONE I=201, DATAPACKING=BLOCK" << endl;

	file.precision(10);
	file.setf(ios::fixed, ios::floatfield);
	for (int i = 0; i < result.size; i++)
	{
		file << i * result.dx - 1.0 << " ";
		file << result.records[i].pressure << " ";
		file << result.records[i].density << " ";
		file << result.records[i].velocity << " ";
		file << endl;
	}

	delete[] result.records;

	return 0;
}

void Solve(Discrete_Records &result,Status &s1, Status &s2) {
	Riemann_Solver solver(s1,s2);
	solver.Solver();
	solver.Get_Result(result);
}