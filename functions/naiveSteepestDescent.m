function [x, fun, dfun, xk] = naiveSteepestDescent(f, F, x0, opts)
% naiveSteepestDescent  Steepest-descent with fixed step or Armijo backtracking
%
%   [x, fun, dfun, xk] = naiveSteepestDescent(f, F, x0, opts)
%
%   Inputs
%     f    : function handle, objective f(x) -> scalar
%     F    : function handle, gradient  F(x) -> column vector
%     x0   : initial guess (column or row vector)
%     opts : (optional) struct with fields:
%            .strategy     'fixed' | 'backtracking'      (default 'backtracking')
%            .stepsize     fixed step size (for 'fixed') (default 1e-2)
%            .alpha        Armijo slope param in (0,0.5] (default 0.2)
%            .beta         backtrack factor in (0,1)     (default 0.5)
%            .maxBacktrack max backtrack inner steps      (default 25)
%            .tol          gradient-norm tolerance        (default 1e-6)
%            .maxIter      maximum iterations             (default 200)
%            .verbose      print per-iteration info       (default true)
%
%   Outputs
%     x    : final iterate
%     fun  : final objective value f(x)
%     dfun : final gradient F(x)
%     xk   : iterate history (each column is an iterate)

    % ---- defaults
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'strategy'),     opts.strategy = 'backtracking'; end
    if ~isfield(opts,'stepsize'),     opts.stepsize = 1e-2; end
    if ~isfield(opts,'alpha'),        opts.alpha = 0.2; end
    if ~isfield(opts,'beta'),         opts.beta = 0.5; end
    if ~isfield(opts,'maxBacktrack'), opts.maxBacktrack = 25; end
    if ~isfield(opts,'tol'),          opts.tol = 1e-6; end
    if ~isfield(opts,'maxIter'),      opts.maxIter = 200; end
    if ~isfield(opts,'verbose'),      opts.verbose = true; end

    % ---- initialization
    x        = x0(:);
    fun      = f(x);
    dfun     = F(x);
    k        = 0;

    n        = numel(x);
    xk       = nan(n, opts.maxIter + 1);
    xk(:,1)  = x;

    if opts.verbose
        fprintf('Iter        f(x)           ||grad||\n');
        fprintf('%4d   %14.8e   %14.8e\n', 0, fun, norm(dfun));
    end

    % ---- main loop
    while (norm(dfun) > opts.tol) && (k < opts.maxIter)
        k = k + 1;

        % Steepest-descent direction
        p = -dfun;

        % Step size selection
        switch lower(opts.strategy)
            case 'fixed'
                t = opts.stepsize;

            case 'backtracking'
                % Armijo: f(x + t p) <= f(x) + alpha*t*grad'*p
                t  = 1.0;
                bt = 0;
                while f(x + t*p) > fun + opts.alpha * t * (dfun.' * p)
                    t  = opts.beta * t;
                    bt = bt + 1;
                    if bt >= opts.maxBacktrack
                        break
                    end
                end

            otherwise
                error('naiveSteepestDescent:BadStrategy', ...
                      'opts.strategy must be ''fixed'' or ''backtracking''.');
        end

        % Update
        x       = x + t*p;
        fun     = f(x);
        dfun    = F(x);
        xk(:,k+1) = x;

        if opts.verbose
            fprintf('%4d   %14.8e   %14.8e\n', k, fun, norm(dfun));
        end
    end

    % trim history
    xk = xk(:,1:k+1);
end
