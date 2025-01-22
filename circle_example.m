clear; close all
%% Initialization

x0.step_length = 0.1;       
x0.xg_sto = {};             
x0.removed = [];           
x0.xg_sto_removed = {}; 

a = 1;
theta = deg2rad(1);
x0.point = [a;0];
x0.point(:,2) = [a*cos(theta);a*sin(theta)];


%% Calling arclength function

func = @(x) circle_eqn(x,a);   
z = arclength8(func,x0,2);  

%% stability analysis
% A_stable=[]; mu_stable = [];
% A_unstable=[]; mu_unstable = [];
% 
% A_arr = z.point(1,:);
% mu_arr = z.point(2,:);
% 
% J_eig = jacobian_eqn_st(func,size(z.point,1));     % calculate jacobian for stability analysis.
% for st_iter = 1:length(z.point)
%     eigen = eig(feval(J_eig,z.point(:,st_iter)));
%     if any(real(eigen)>0)
%         A_unstable = [A_unstable,A_arr(st_iter)];
%         mu_unstable = [mu_unstable,mu_arr(st_iter)];
%     else
%         A_stable = [A_stable,A_arr(st_iter)];
%         mu_stable = [mu_stable,mu_arr(st_iter)];
%     end
% end

%% plotting

% figure(2)
% scatter(mu_stable,A_stable,"filled","blue","DisplayName","Stable")
% hold on
% scatter(mu_unstable,A_unstable,"filled","red","DisplayName","Unstable")
% xlabel("$\mu_{0}$","FontSize",20,"Interpreter","latex")
% ylabel("$A$","FontSize",20,"Interpreter","latex")
% legend("Location","best","Interpreter","latex","FontSize",18)
% ax = gca;
% ax.FontSize = 18;


%% plotting

figure(2)
scatter(z.point(1,:),z.point(2,:),"filled","blue")
xlabel("$x$","FontSize",20,"Interpreter","latex")
ylabel("$y$","FontSize",20,"Interpreter","latex")
ax = gca;
ax.FontSize = 18;


%% Functions
function z = circle_eqn(x,a)

z = x(1)^2 + x(2)^2 - a^2;

end