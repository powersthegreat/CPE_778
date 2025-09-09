clc; close;

addpath('../../functions');

% ----- define objective, gradient, Hessian
f = @(v) (1 - v(1)).^2 + 100*(v(2) - v(1).^2).^2;

F = @(v) [ -2*(1 - v(1)) - 400*v(1)*(v(2) - v(1).^2) ; ...
            200*(v(2) - v(1).^2) ];

H = @(v) [ 2 - 400*(v(2) - 3*v(1).^2),   -400*v(1); ...
          -400*v(1),                      200      ];

% ----- initial guess
x0 = [-1.2; 1.0];

% ----- run Newton
opts = struct('alpha',0.2,'beta',0.5,'tol',1e-8,'maxiter',200,'verbose',true);
[x, info] = amijoNewton(f, F, H, x0, opts);

fprintf('\nResult:\n');
fprintf('x* = [% .8f, % .8f]^T\n', x(1), x(2));
fprintf('f(x*) = %.3e,  ||grad|| = %.3e,  iters = %d, stop=%s\n', ...
        f(x), norm(F(x)), info.iters, info.stops);

% ----- simple checks
assert(norm(F(x)) < 1e-6, 'Gradient tolerance not met.');
assert(norm(x - [1;1]) < 1e-3, 'Did not converge to [1;1] within tolerance.');

% ----- plot contour and path
plot_rosenbrock_with_path(f, info.xs);


function plot_rosenbrock_with_path(f, xs)
    % grid for contour
    [X, Y] = meshgrid(linspace(-2, 2, 200), linspace(-1, 3, 200));
    Z = arrayfun(@(i,j) f([i;j]), X, Y);

    figure; hold on; grid on;
    contour(X, Y, Z, 10.^(-2:2), 'ShowText','on'); % log-scaled levels
    plot(xs(1,:), xs(2,:), 'o-','LineWidth',1.5);
    plot(1,1,'kp','MarkerSize',12,'MarkerFaceColor','y'); % optimum
    xlabel('x'); ylabel('y'); title('myNewton on Rosenbrock');
    legend('Contours','Iterates','[1,1]','Location','best');
end
