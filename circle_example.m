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
z = arclength8(func,x0,100);  

%% plotting

figure(2)
scatter(z.point(1,:),z.point(2,:),"filled","blue")
xlabel("$x$","FontSize",20,"Interpreter","latex")
ylabel("$y$","FontSize",20,"Interpreter","latex")
ax = gca;
ax.XLim = [-1 1];
ax.YLim = [-1 1];
ax.FontSize = 18;
axis equal

%% Functions
function z = circle_eqn(x,a)

z = x(1)^2 + x(2)^2 - a^2;

end