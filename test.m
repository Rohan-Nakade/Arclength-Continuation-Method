clear; close all
% %% Initialization
% % create a struct with required parameters as fields of the struct.
% 
% x0.step_length = 0.1;       % initial step length
% x0.xg_sto = {};             % cell to store guess values (can be used for debugging)
% x0.removed = [];            % array which will contain the point that is removed as per the algorithm
% x0.xg_sto_removed = {};     % guess values removed corresponding to above point
% 
% % finding the first point from where the continuation will start using
% % ode45 (can also use RK4).
% % Here we will start from mu = 2 and go backward. This is because when A
% % becomes zero, the jacobians will become non-invertible and issues will
% % occur during continuation.
% 
% xinit = 0.1;
% tspan = 0:0.1:100;
% indx = ceil(2*length(tspan)/3);
% 
% mu = 2;
% [t,x_final] = ode45(@(t,x) vanderpol(x,mu),tspan,xinit);
% x_rms = rms(x_final(indx:end,1));
% x0.point(1,1) = x_rms;          % 1st row will contain the value of A
% x0.point(2,1) = mu;             % 2nd row will contain the value of control parameter (mu)
%                                 % If have more variables, then the control parameter should be in last row
% 
% mu = 1.99;
% [t,x_final] = ode45(@(t,x) vanderpol(x,mu),tspan,xinit);
% x_rms = rms(x_final(indx:end,1));
% x0.point(1,2) = x_rms;
% x0.point(2,2) = mu;
% 
% % Here we are calling ode45 two times to specify the direction in which we
% % want to continue i.e. from mu=2 to mu=1.99 (this will go in decreasing mu direction)
%% Calling arclength function

func = @(x) vanderpol(x,x(2));   % create a new function with control parameter as a variable.
% z = arclength8(func,x0,100);     % 100 is the number of loops. The arclength function will give that many points.

% data type of z will be struct which will have same fields as x0.

%% stability analysis

J_eig = @(x) jacobian_eqn_st(func,2);     % calculate jacobian for stability analysis.
eigen = @(x) eig(feval(J_eig,x));
    


%% Functions
function z = vanderpol(x,mu0)
% mu6 = 0.01; mu4 = -0.6; mu2 = 7; 
mu6 = 0; mu4 = 3; mu2 = -10; 
omega = 220*2*pi; gamma = 0;
A = x(1);

Adot = mu0*A/2 - mu2*A^3/8 - mu4*A^5/16 - 5*mu6*A^7/128 + gamma/(4*omega^2*A);

z = Adot;
end