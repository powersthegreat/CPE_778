clc; clear;

% Add path for custom functions (if needed)
addpath('../../functions');

% Case 1: 2-var quadratic
syms x y
f  = x^2 + x*y + 3*y;
g_expected = [2*x + y; x + 3];
g = computeGradient(f, [x y]);
assert(isequal(simplify(g - g_expected), sym(zeros(2,1))), 'Case 1 failed.');
fprintf('  Case 1: PASS\n');

% Numeric spot-check
gv = subs(g, [x y], [2 1]);
assert(isequal(simplify(gv - [5; 5]), sym(zeros(2,1))), 'Case 1 numeric failed.');

% Case 2: 3-var mixed trig/poly
syms x y z
f2 = sin(x*y) + z^3 + x*z;
g2_expected = [ y*cos(x*y) + z;  x*cos(x*y);  3*z^2 + x ];
g2 = computeGradient(f2, [x y z]);
assert(isequal(simplify(g2 - g2_expected), sym(zeros(3,1))), 'Case 2 failed.');
fprintf('  Case 2: PASS\n');

% Numeric spot-check
gv2 = subs(g2, [x y z], [1 2 0]);
assert(isequal(simplify(gv2 - [2*cos(2) + 0; 1*cos(2); 0 + 1]), sym(zeros(3,1))), ...
       'Case 2 numeric failed.');

% Case 3: Default variable ordering (symvar) consistency
% Note: symvar order may differ across versions/files; we only check that
% calling computeGradient without vars matches gradient(f, symvar(f)).
syms a b c
f3 = a^2 + 2*b + 3*c + a*b;
g3 = computeGradient(f3);                 % our function
g3_direct = gradient(f3, symvar(f3));     % MATLAB direct
assert(isequal(simplify(g3 - g3_direct), sym(zeros(size(g3)))), ...
       'Case 3 symvar ordering mismatch.');
fprintf('  Case 3: PASS\n');

fprintf('All computeGradient tests PASSED.\n');
