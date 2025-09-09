Variables
    x1, x2, z;

Equations
    obj;

* Objective function (quadratic)
obj.. z =e= 5*sqr(x1) + sqr(x2) + 2*x1*x2 - 12*x1 - 4*x2 + 8;

Model quad /all/;

* Starting point 1
x1.l = 0;
x2.l = -2;

Solve quad using nlp minimizing z;

display x1.l, x2.l, z.l;

Parameter
    results(*);

results('x1')        = x1.l;
results('x2')        = x2.l;
results('z')         = z.l;
results('ModelStat') = quad.modelstat;
results('SolveStat') = quad.solvestat;
results('CPUtime')   = timeelapsed;

display results;


* Starting point 2
x1.l = 0;
x2.l = -4;

Solve quad using nlp minimizing z;

display x1.l, x2.l, z.l;

Parameter
    results(*);

results('x1')        = x1.l;
results('x2')        = x2.l;
results('z')         = z.l;
results('ModelStat') = quad.modelstat;
results('SolveStat') = quad.solvestat;
results('CPUtime')   = timeelapsed;

display results;

* Starting point 3
x1.l = 2;
x2.l = -2;

Solve quad using nlp minimizing z;

display x1.l, x2.l, z.l;

Parameter
    results(*);

results('x1')        = x1.l;
results('x2')        = x2.l;
results('z')         = z.l;
results('ModelStat') = quad.modelstat;
results('SolveStat') = quad.solvestat;
results('CPUtime')   = timeelapsed;

display results;