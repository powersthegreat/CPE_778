Set
   i / 1, 2, 3, 4 /;

Positive Variable x(i);
Variable z;

Equation
   cost
   eq1
   eq2
   eq3;

* Objective function
cost.. z =e= 8*x('1') + 9*x('2') + 5*x('3');

* Constraints
eq1.. x('1') + x('2') + x('3') =l= 2;
eq2.. 2*x('1') + 3*x('2') + 4*x('3') =l= 3;
eq3.. 6*x('1') + 6*x('2') + 2*x('3') =l= 8;

Model test / all /;

solve test using lp minimizing z;

display x.l, x.m;


$ontext

Solution Output:

Tried aggregator 1 time.
LP Presolve eliminated 3 rows and 4 columns.
All rows and columns eliminated.
Presolve time = 0.00 sec. (0.00 ticks)

--- LP status (1): optimal.
--- Cplex Time: 0.00sec (det. 0.00 ticks)


Optimal solution found
Objective:            0.000000

$offtext
