clc; clear;

% Add path for custom functions (if needed)
addpath('../../functions');

% Define symbolic variables
syms d n beta;

% Define function
f = (beta/d + d^n)^n;

% Example: calculate 1st derivative
F1 = nthDerivative(f, d, 1);
disp('First derivative:');
disp(F1);

% Example: calculate 2nd derivative
F2 = nthDerivative(f, d, 2);
disp('Second derivative:');
disp(F2);
