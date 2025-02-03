% This function calculates the solution point after the guess is updated
% from the arclength8.m. This function uses the Newton Raphson algorithm to
% get convergence and updates the guesses on a circle of radius equal to
% the step length.

% This function takes the inputs:
%   ds - step length.
%   xstore - the array of calculated solution points.
%   f - function for which solution branch needs to be calculated.
%   Jx - Jacobian of function f.
%   xg - the updated guess from the arclength8.m file.
%   a - relaxation parameter for ensuring smooth convergence.

% The output contains:
%   xg - The updated guess from NR method.
%   successfull - tells if the convergence was met or not (1 if convergence and 0 otherwise)
%   xg_store - array of all the guesses obtained in the iterations.
%   A - used for debugging

function [xg,successful,xg_store, A] = cont4(ds,xstore,Jx,f,xg,a)
    xg_store = xg;
    xn = xstore.point(:,end);
    dx = xg - xn;

    tol = norm(dx)*1e-2;    % tolerance for convergence.
    maxiter = 2e3;          % max number of iterations.
    numiter = 0;

    while norm(dx) > tol && numiter < maxiter
        J = feval(Jx,xg);
        A1 = zeros(1,length(J(1,:)));
        for jj = 1:length(J(1,:))
            B1 = J(end,:);
            B1(jj) = 1;
            A1(jj) = prod(B1);
        end
        A = [J ; A1];
        B = [-feval(f,xg) ; 0];
        dx = A\B;
        xg = xg + a*dx;

        % force the solution on the circle.
        Ds = norm(xg - xn);
        Dx = xg - xn;
        xg = xn + Dx*ds/Ds;
        
        numiter = numiter + 1;
        xg_store = [xg_store,xg];

    end

    % check if the convergence was successfull or not.
    if norm(dx)<=tol
        successful = 1;
    else
        successful = 0;
    end
    
end