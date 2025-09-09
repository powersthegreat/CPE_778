
Set
   i / 1, 2 /
   j / 1, 2, 3, 4 /;
          
Table c(i,j) 
              1        2       3       4
   1         2.5      1.7     1.8     2.0
   2         2.5      1.8     1.4     2.2;

Variable x(i,j);

Positive Variable x;

Equation
   cost
   eq1(j)
   eq2(j)
   eq3(i)
   eq4(i)
   eq5(i)
   eq6(j) ;

cost..      z =e= sum((i,j), c(i,j)*x(i,j));

eq1(j).. sum(j, x(1,j)) =l= 150;
eq2(j).. sum(j, x(2,j)) =l= 200;

eq3(i).. sum(i, x(i, 1)) =g= 100;
eq4(i).. sum(i, x(i, 2)) =g= 80;
eq5(i).. sum(i, x(i, 3)) =g= 60;
eq6(i).. sum(i, x(i, 4)) =g= 900;

Model test / all /;

solve test using lp minimizing z;

display x.l, x.m;

