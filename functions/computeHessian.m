function H = computeHessian(f, vars)
% computeHessian  Hessian (matrix of second partials) of a symbolic scalar f
%
%   H = computeHessian(f, vars)
%
%   INPUTS
%     f    : symbolic scalar expression (e.g., f = x^2 + x*y)
%     vars : (optional) symbolic variable vector that orders the Hessian
%            rows/cols, e.g., vars = [x y]. If omitted, symvar(f) is used.
%
%   OUTPUT
%     H    : symbolic Hessian matrix (size length(vars)-by-length(vars))
%
%   EXAMPLES
%     syms x y
%     f = x^2 + x*y + 3*y;
%     H = computeHessian(f, [x y])
%     % returns: [ 2  1
%     %            1  0 ]
%
%   TIP
%     Evaluate at a point with SUBS:
%       subs(H, [x y], [2 1])

    arguments
        f (1,1) sym
        vars sym = symvar(f)
    end

    if ~isvector(vars) || ~all(arrayfun(@(s) isa(s,'sym'), vars))
        error('computeHessian:VarsMustBeSymVector', ...
              'vars must be a vector of symbolic variables, e.g., [x y z].');
    end

    H = hessian(f, vars);
end
