function [x, info] = amijoNewton(f, F, H, x0, opts)
% amijoNewton  Newton's method with Armijo backtracking line search
%
%   [x, info] = myNewton(f, F, H, x0, opts)
%
%   Inputs
%     f    : handle, objective f(x) -> scalar
%     F    : handle, gradient  F(x) -> column vector
%     H    : handle, Hessian   H(x) -> matrix
%     x0   : column vector, initial guess
%     opts : (optional) struct with fields:
%            .alpha     Armijo slope parameter in (0, 0.5],   default 0.2
%            .beta      backtracking factor in (0,1),         default 0.5
%            .tol       gradient norm tolerance,              default 1e-6
%            .maxiter   maximum iterations,                   default 100
%            .verbose   print per-iter info (true/false),     default true
%
%   Outputs
%     x    : final iterate
%     info : struct with fields:
%            .iters      number of iterations
%            .grad_norms gradient norms per iteration
%            .fvals      objective values per iteration
%            .xs         iterates (each column is x_k)
%            .stops      string stop reason

    % ---- options & defaults
    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'alpha'),   opts.alpha   = 0.2;  end   % "a" in the excerpt
    if ~isfield(opts,'beta'),    opts.beta    = 0.5;  end   % "beta" in the excerpt
    if ~isfield(opts,'tol'),     opts.tol     = 1e-6; end
    if ~isfield(opts,'maxiter'), opts.maxiter = 100;  end
    if ~isfield(opts,'verbose'), opts.verbose = true; end

    alpha   = opts.alpha;
    beta    = opts.beta;
    tol     = opts.tol;
    maxiter = opts.maxiter;

    % ---- initialization (STEP 1 in the excerpt)
    x    = x0(:);
    fun  = f(x);
    grad = F(x);
    Hk   = H(x);

    % Newton step (STEP 2 in the excerpt)
    p = compute_newton_step(Hk, grad);

    % bookkeeping
    info.iters      = 0;
    info.grad_norms = norm(grad);
    info.fvals      = fun;
    info.xs         = x;

    if opts.verbose
        fprintf('Iter        f(x)           ||grad||\n');
        fprintf('%4d   %14.8e   %14.8e\n', 0, fun, norm(grad));
    end

    % ---- main loop (STEP 3 in the excerpt)
    while (norm(grad) > tol) && (info.iters < maxiter)
        info.iters = info.iters + 1;

        % Armijo backtracking: find n s.t.
        % f(x + beta^n p) <= f(x) + alpha * beta^n * grad' * p
        n = 0;
        t = beta^n;
        while f(x + t*p) > fun + alpha * t * (grad.'*p)
            n = n + 1;
            t = beta^n;
            if n == 50          % hard cap to avoid infinite loop
                break;
            end
        end

        % Update iterate (STEP 2 again)
        x   = x + t*p;
        fun = f(x);
        grad = F(x);
        Hk   = H(x);
        p    = compute_newton_step(Hk, grad);

        % log
        info.grad_norms(end+1) = norm(grad);
        info.fvals(end+1)      = fun;
        info.xs(:, end+1)      = x;

        if opts.verbose
            fprintf('%4d   %14.8e   %14.8e\n', info.iters, fun, norm(grad));
        end
    end

    if norm(grad) <= tol
        info.stops = 'gradient_tol_met';
    else
        info.stops = 'maxiter_reached';
    end
end

% ---- helper: robust Newton step (fallback to -grad if Hessian is bad)
function p = compute_newton_step(Hk, g)
    % Try Newton direction
    % If Hessian is singular or does not give a descent direction,
    % fall back to steepest descent.
    try
        % symmetric solve is more stable when Hk is SPD
        p = - (Hk \ g);
        if g.'*p >= 0         % not a descent direction
            p = -g;
        end
    catch
        p = -g;               % singular/ill-conditioned Hessian
    end
end
