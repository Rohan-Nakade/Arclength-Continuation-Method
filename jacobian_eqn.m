function z = jacobian_eqn(f,n)
x = sym("x",[n,1]);
expression = feval(f,x);
jac = jacobian(expression,x);

z = matlabFunction(jac,"Vars",{x});
end