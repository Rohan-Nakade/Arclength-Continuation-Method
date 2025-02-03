% This function gives the Jacobian of f for computing stability.
% J(X) = d(f(X))/dX
% x = [X;K] - has dimensionality of n and X has dimensionality of n-1.

function z = jacobian_eqn_st(f,n)
x = sym("x",[n,1]);
expression = feval(f,x);
jac = jacobian(expression,x(1:n-1));

z = matlabFunction(jac,"Vars",{x});
end