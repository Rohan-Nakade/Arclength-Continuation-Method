% This function generates a function for the Jacobian equation for f.
% This function takes inputs:
%   f - function for which solution branch is to be computed
%   n - dimensionality 'x' which contains the dependable variables X and
%       the control/bifurcation parameter

% x = [X;K] - has dimensionality of n
% f = f(X) - has dimensionality of n-1
% J = df/dx

function z = jacobian_eqn(f,n)
x = sym("x",[n,1]);
expression = feval(f,x);
jac = jacobian(expression,x);

z = matlabFunction(jac,"Vars",{x});
end