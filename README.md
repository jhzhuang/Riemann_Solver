# Riemann_Solver
  This solver deals with the one-dimensional Riemann problem in gas dynamics with one discontinuity at x = 0 
and gives the result of analytical solution at disctrete nodes. In my program, there is no necessarity 
to seperate the initial discontinuity problem into five cases. The solver decide whether a wave is a shock 
wave or expansion wave by simply comparing the pressure before it and after it, and this makes it possible to automaticly 
compute parameters in certain zones using respective formulas for shock wave or expansion wave. The file example.cpp 
shows how to use this solver.

