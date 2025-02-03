% This script gives the solution branches for the Vanderpol oscillator.

clear; close all
%% Initialization
% create a struct with required parameters as fields of the struct.

x0.step_length = 0.1;       % initial step length
x0.xg_sto = {};             % cell to store guess values (can be used for debugging)
x0.removed = [];            % array which will contain the point that is removed as per the algorithm
x0.xg_sto_removed = {};     % guess values removed corresponding to above point

% finding the first point from where the continuation will start using
% ode45 (can also use RK4).
% Here we will start from mu = 2 and go backward. This is because when A
% becomes zero, the jacobians will become non-invertible and issues will
% occur during continuation.

xinit = 0.1;
tspan = 0:0.1:100;
indx = ceil(2*length(tspan)/3);

mu = 2;
[t,x_final] = ode45(@(t,x) vanderpol(x,mu),tspan,xinit);
x_rms = rms(x_final(indx:end,1));
x0.point(1,1) = x_rms;          % 1st row will contain the value of A
x0.point(2,1) = mu;             % 2nd row will contain the value of control parameter (mu)
                                % If have more variables, then the control parameter should be in last row

mu = 1.99;
[t,x_final] = ode45(@(t,x) vanderpol(x,mu),tspan,xinit);
x_rms = rms(x_final(indx:end,1));
x0.point(1,2) = x_rms;
x0.point(2,2) = mu;

% Here we are calling ode45 two times to specify the direction in which we
% want to continue i.e. from mu=2 to mu=1.99 (this will go in decreasing mu direction)
%% Calling arclength function

func = @(x) vanderpol(x,x(2));   % create a new function with control parameter as a variable.
z = arclength8(func,x0,150);     % 100 is the number of loops. The arclength function will give that many points.

% data type of z will be struct which will have same fields as x0.

%% stability analysis
A_stable=[]; mu_stable = [];
A_unstable=[]; mu_unstable = [];

A_arr = z.point(1,:);
mu_arr = z.point(2,:);

J_eig = jacobian_eqn_st(func,size(z.point,1));     % calculate jacobian for stability analysis.
for st_iter = 1:length(z.point)
    eigen = eig(feval(J_eig,z.point(:,st_iter)));
    if any(real(eigen)>0)
        A_unstable = [A_unstable,A_arr(st_iter)];
        mu_unstable = [mu_unstable,mu_arr(st_iter)];
    else
        A_stable = [A_stable,A_arr(st_iter)];
        mu_stable = [mu_stable,mu_arr(st_iter)];
    end
end

%% plotting

figure(2)
scatter(mu_stable,A_stable,"filled","blue","DisplayName","Stable")
hold on
scatter(mu_unstable,A_unstable,"filled","red","DisplayName","Unstable")
xlabel("$\mu_{0}$","FontSize",20,"Interpreter","latex")
ylabel("$A$","FontSize",20,"Interpreter","latex")
legend("Location","best","Interpreter","latex","FontSize",18)
ax = gca;
ax.FontSize = 18;


%% Functions
function z = vanderpol(x,mu0)
% mu6 = 0.01; mu4 = -0.6; mu2 = 7; 
mu6 = 0; mu4 = 3; mu2 = -10; 
omega = 220*2*pi; gamma = 0;
A = x(1);

Adot = mu0*A/2 - mu2*A^3/8 - mu4*A^5/16 - 5*mu6*A^7/128 + gamma/(4*omega^2*A);

z = Adot;
end