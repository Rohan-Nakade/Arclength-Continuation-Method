function z = jacobian_eqn_st(f,n)
x = sym("x",[n,1]);
expression = feval(f,x);
jac = jacobian(expression,x(1:n-1));

z = matlabFunction(jac,"Vars",{x});
end