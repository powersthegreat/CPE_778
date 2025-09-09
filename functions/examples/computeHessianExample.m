clc; clear;

% Add path for custom functions (if needed)
addpath('../../functions');


% Case 1: 2-var quadratic
syms x y
f  = x^2 + x*y + 3*y;
H_expected = [2 1; 1 0];
H = computeHessian(f, [x y]);
assert(isequal(simplify(H - H_expected), sym(zeros(2))), 'Case 1 failed.');
fprintf('  Case 1: PASS\n');

% Case 2: 3-var mixed trig/poly
syms x y z
f2 = sin(x*y) + z^3 + x*z;
% d/dx:  y*cos(x*y) + z
% d/dy:  x*cos(x*y)
% d/dz:  3*z^2 + x
% Hessian entries:
% Hxx = -y^2*sin(x*y), Hxy = cos(x*y) - x*y*sin(x*y), Hxz = 1
% Hyx = same as Hxy,   Hyy = -x^2*sin(x*y),           Hyz = 0
% Hzx = 1,             Hzy = 0,                       Hzz = 6*z
H2_expected = [ -y^2*sin(x*y),   cos(x*y) - x*y*sin(x*y), 1;
                 cos(x*y) - x*y*sin(x*y), -x^2*sin(x*y),  0;
                 1,                0,                     6*z ];

H2 = computeHessian(f2, [x y z]);
assert(isequal(simplify(H2 - H2_expected), sym(zeros(3))), 'Case 2 failed.');
fprintf('  Case 2: PASS\n');

% Numeric spot-check
H2v = subs(H2, [x y z], [1 2 3]);
H2v_exp = subs(H2_expected, [x y z], [1 2 3]);
assert(isequal(simplify(H2v - H2v_exp), sym(zeros(3))), 'Case 2 numeric failed.');

% Case 3: Default variable ordering (symvar) consistency
syms a b c
f3 = a^2 + 2*b + 3*c + a*b;
H3 = computeHessian(f3);
H3_direct = hessian(f3, symvar(f3));
assert(isequal(simplify(H3 - H3_direct), sym(zeros(size(H3)))), ...
       'Case 3 symvar ordering mismatch.');
fprintf('  Case 3: PASS\n');

fprintf('All computeHessian tests PASSED.\n');

