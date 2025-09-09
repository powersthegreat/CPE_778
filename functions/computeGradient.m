function g = computeGradient(f, vars)
% computeGradient  Gradient of a multivariate symbolic function
%
%   g = computeGradient(f, vars)
%
%   INPUTS
%     f    : symbolic scalar expression (e.g., f = x^2 + x*y)
%     vars : (optional) symbolic variable vector that orders the gradient,
%            e.g., vars = [x y]. If omitted, symvar(f) is used (MATLAB’s
%            default variable ordering), which may differ from your intent.
%
%   OUTPUT
%     g    : column vector of partial derivatives (symbolic)
%
%   EXAMPLES
%     syms x y
%     f = x^2 + x*y + 3*y;
%     g = computeGradient(f, [x y])
%     % returns: [ 2*x + y;  x + 3 ]
%
%   TIP
%     To evaluate at a point, use SUBS:
%       subs(g, [x y], [2 1])   % -> [2*2 + 1;  2 + 3] = [5; 5]

    arguments
        f (1,1) sym
        vars sym = symvar(f)  % default: MATLAB’s inferred variable order
    end

    % Ensure 'vars' is a vector of symbols
    if ~isvector(vars) || ~all(arrayfun(@(s) isa(s,'sym'), vars))
        error('computeGradient:VarsMustBeSymVector', ...
              'vars must be a vector of symbolic variables, e.g., [x y z].');
    end

    % Compute gradient (column vector); gradient(f, vars) respects ordering
    g = gradient(f, vars);
end
