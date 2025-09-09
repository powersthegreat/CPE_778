function [x, fun, dfun, xk] = naiveNewton(f, F, H, x0)
% NaiveNewton  Newton’s method with backtracking line search
%
%   [x, fun, dfun, xk] = NaiveNewton(f, F, H, x0)
%
%   Inputs:
%     f  : function handle (objective)
%     F  : function handle (first derivative)
%     H  : function handle (second derivative)
%     x0 : initial guess
%
%   Outputs:
%     x   : final solution
%     fun : objective value at solution
%     dfun: derivative value at solution
%     xk  : iterates history

    % Parameters
    maxiter       = 100;   % maximum iterations
    maxbacktrack  = 10;    % maximum backtracking steps
    tol           = 1e-6;  % stopping tolerance

    % Initialization
    k  = 0;
    x  = x0;
    xk = x0;

    fprintf('Iter         x             f(x)           ||g(x)||\n');

    while norm(feval(F, x)) > tol && k < maxiter
        % Evaluate function and derivatives
        fun  = feval(f, x);
        dfun = feval(F, x);
        ddfun = feval(H, x);

        % Newton step
        FullNewtonStep = -dfun / ddfun;

        % Backtracking line search
        n     = 0;
        newx  = x + (1/2)^n * FullNewtonStep;

        while feval(f, newx) > fun && n < maxbacktrack
            n    = n + 1;
            newx = x + (1/2)^n * FullNewtonStep;
            fprintf('  Backtracking with n = %4u\n', n);
        end

        if n == maxbacktrack
            k = maxiter;   % terminate if max backtracking exceeded
        end

        % Update
        x  = newx;
        k  = k + 1;
        xk(:, k+1) = x;

        % Print iteration info
        fprintf('%4u   %14.8f   %14.8f   %14.8f\n', ...
                k, x, fun, norm(dfun));
    end

    % Final check
    if ddfun > eps
        fprintf('Possible local minimum\n');
    else
        fprintf('Newton has failed\n');
    end
end
