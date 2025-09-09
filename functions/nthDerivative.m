function Fn = nthDerivative(f, x, n)
% nthDerivative Compute the nth derivative of a symbolic function
%
%   Fn = nthDerivative(f, x, n) computes the nth derivative of the
%   symbolic expression f with respect to the symbolic variable x.
%
%   Inputs:
%     f : symbolic expression (function of x)
%     x : symbolic variable to differentiate with respect to
%     n : order of derivative (integer)
%
%   Output:
%     Fn : symbolic expression of the nth derivative

    arguments
        f sym
        x sym
        n (1,1) {mustBeInteger, mustBeNonnegative}
    end

    Fn = simplify(diff(f, x, n));
end