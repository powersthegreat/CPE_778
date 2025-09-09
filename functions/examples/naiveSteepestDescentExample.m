clc; close;

addpath('../../functions');

% -------- Example 1: Quadratic, FIXED stepsize  ------------------------
% f(x) = 1/2 x^T Q x + b^T x  (Q SPD)
Q = diag([1, 10]);
b = [-1; 2];

f_quad = @(x) 0.5 * (x.' * (Q*x)) + b.'*x;
F_quad = @(x) Q*x + b;

x0 = [5; -3];
L  = max(eig(Q));       % Lipschitz constant of grad
eta = 1/L;              % optimal fixed step for gradient descent on quadratics

opts1 = struct('strategy','fixed','stepsize',eta, ...
               'tol',1e-10,'maxIter',500,'verbose',true);

[x1, f1, g1, xk1] = naiveSteepestDescent(f_quad, F_quad, x0, opts1);

fprintf('\nQuadratic w/ fixed step:\n');
fprintf('  x*   = [% .6f, % .6f]^T\n', x1(1), x1(2));
fprintf('  f(x*)= %.3e,  ||grad||= %.3e\n', f1, norm(g1));

% closed-form solution for comparison: Qx + b = 0 -> x* = -Q^{-1} b
x_star = -Q \ b;
assert(norm(x1 - x_star) < 1e-6, 'Fixed-step did not reach the quadratic optimum.');
assert(norm(g1) < 1e-8, 'Gradient not sufficiently small for quadratic.');

% -------- Example 2: Rosenbrock, BACKTRACKING --------------------------
% f(x,y) = (1-x)^2 + 100 (y - x^2)^2, min at [1;1]
f_ros = @(v) (1 - v(1)).^2 + 100*(v(2) - v(1).^2).^2;
F_ros = @(v) [ -2*(1 - v(1)) - 400*v(1)*(v(2)-v(1).^2) ; ...
                200*(v(2) - v(1).^2) ];

x0 = [-1.2; 1.0];

opts2 = struct('strategy','backtracking','alpha',0.2,'beta',0.5, ...
               'tol',1e-6,'maxIter',5000,'verbose',true);

[x2, f2, g2, xk2] = naiveSteepestDescent(f_ros, F_ros, x0, opts2);

fprintf('\nRosenbrock w/ backtracking:\n');
fprintf('  x*   = [% .6f, % .6f]^T\n', x2(1), x2(2));
fprintf('  f(x*)= %.3e,  ||grad||= %.3e\n', f2, norm(g2));

assert(norm(g2) < 1e-4, 'Backtracking descent did not reduce gradient enough.');
assert(norm(x2 - [1;1]) < 1e-2, 'Did not get sufficiently close to [1;1].');

% Visualize path (optional)
plot_rosenbrock_path(f_ros, xk2);


function plot_rosenbrock_path(f, xs)
    [X, Y] = meshgrid(linspace(-2,2,200), linspace(-1,3,200));
    Z = arrayfun(@(i,j) f([i;j]), X, Y);
    figure; hold on; grid on;
    contour(X, Y, Z, 10.^(-2:2), 'ShowText','on');
    plot(xs(1,:), xs(2,:), 'o-','LineWidth',1.5);
    plot(1,1,'kp','MarkerSize',12,'MarkerFaceColor','y');
    xlabel('x'); ylabel('y'); title('Steepest Descent (backtracking) on Rosenbrock');
    legend('Contours','Iterates','[1,1]','Location','best');
end
