% This function calculates the solution branches for the function 'f'.
% x0 a struct which contains the initial condition for the following fields:
%       x0.step_length - contains the initial steplength used for arclength continuation
%       x0.xg_sto - cell to store guess values (can be used for debugging)
%       x0.removed - array which will contain the point that is removed as per the algorithm
%       x0.xg_sto_removed - guess values removed corresponding to above point

%   The algorithm used can be found in DDE-BIFTOOL repository on GitHub: https://github.com/DDE-BifTool/DDE-Biftool

function z = arclength8(f,x0,num)

growth_factor = 1.2;        % after each successful convergence, the steplength is increased by this factor.
relaxation_param = 0.1;     % used for smooth convergence and remove the chance of oscillation about the solution point.
relaxation_param_new = 0.1;
step_min = 0.001;           % min value of step length
step_max = 0.1;             % max value of step length
xstore = x0;

Jx = jacobian_eqn(f,length(x0.point));  % generates a function for the Jacobian of 'f'.

for ii = 1:num
    disp(ii)
    ds = max(step_min,xstore.step_length(end));
    if ii == 1
        xg = xstore.point(:,end) + ds*secant(xstore.point);   % first guess
        [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param); % solution using Newton-Raphson method
        v = secant(xstore.point(:,end),xg_new);
        xstore.secant = v;
        xstore.point(:,end) = [];
    else
        slope = secant(xstore.point);
        xg = xstore.point(:,end) + ds*slope;
        [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param);
        v = secant(xstore.point(:,end),xg_new);
        while ~success && (ds>step_min)     % reduce step length until convergence is met or steplength is smaller than min value of step length
            ds = ds/2;
            relaxation_param_new = relaxation_param_new*0.9;
            xg = xstore.point(:,end) + ds*slope;
            [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param_new);
            v = secant(xstore.point(:,end),xg_new);
        end
        relaxation_param_new = relaxation_param;
    end

    V = xstore.secant(:,end);           % store the slope information at the calculated solution point.

    % if convergence is met, then check the direction of new point with
    % respect to previous one. If direction is forward then add the point.
    % Otherwise, store the point in between the previous-1 and previous.
    if success
        relaxation_param_new = relaxation_param;
        if ds < step_max
            ds = min(ds*growth_factor,step_max);
        end
        if V'*v > 0
            xstore.point = [xstore.point , xg_new];
            xstore.step_length = [xstore.step_length,ds];
            xstore.xg_sto = [xstore.xg_sto, guess];
            xstore.secant = [xstore.secant,v];
        elseif V'*v < 0
            xstore.point(:,end+1) = xstore.point(:,end);
            xstore.point(:,end-1) = xg_new;
            xstore.step_length(:,end+1) = ds*growth_factor;
            xstore.step_length(:,end-1) = ds;
            xstore.xg_sto(end+1) = xstore.xg_sto(end);
            xstore.xg_sto(end-1) = {guess};
            xstore.secant(:,end+1) = xstore.secant(:,end);
            xstore.secant(:,end-1) = secant(xstore.point);
        end
    end

    % If convergence is not met then follow the algorithm given by
    % DDE-BIFTOOL.
    % Calculate a new point in backward direction at half steplength. If no
    % convergence if met at this point, remove the previous point.
    % Then check the direction of this half steplength point and add as per
    % the direction in the array of solution points.
    if ~success
        ds = -ds/2;
        xg = xstore.point(:,end) + ds*slope;
        [xg_new,success,guess] = cont4(-ds,xstore,Jx,f,xg,relaxation_param);
        check_v = secant(xstore.point(:,end),xg_new);

        if ~success && ~isempty(xstore.point)
            xstore.removed = [xstore.removed,xstore.point(:,end)];
            xstore.point(:,end)=[];
            xstore.step_length(end) = [];
            xstore.step_length(end) = 0.01;
            xstore.xg_sto(end) = [];
            xstore.secant(:,end) = [];
        elseif success && norm(xstore.point(end-1) - xstore.point(end-2))>0
            if check_v'*V > 0
                xstore.point = [xstore.point , xg_new];
                xstore.step_length = [xstore.step_length,ds];
                xstore.xg_sto = [xstore.xg_sto, guess];
                xstore.secant = [xstore.secant,check_v];
            else
                xstore.point(:,end+1) = xstore.point(:,end);
                xstore.point(:,end-1) = xg_new;
                xstore.step_length(:,end+1) = -ds*growth_factor;
                xstore.step_length(:,end-1) = -ds;
                xstore.xg_sto(end+1) = xstore.xg_sto(end);
                xstore.xg_sto(end-1) = {guess};
                xstore.secant(:,end+1) = xstore.secant(:,end);
                xstore.secant(:,end-1) = secant(xstore.point);
            end
        end
    end
end

z = xstore;

end