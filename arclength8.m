function [z,test_arr] = arclength8(f,x0,num)
global test_arr
growth_factor = 1.2;
relaxation_param = 0.1;
relaxation_param_new = 0.4;
step_min = 0.001;
step_max = 0.1;
xstore = x0;
Jx = jacobian_eqn(f,length(x0.point));
for ii = 1:num
    disp(ii)
    ds = max(step_min,xstore.step_length(end));
    if ii == 1
        xg = xstore.point(:,end) + ds*secant(xstore.point);   %first guess
        [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param);
        v = secant(xstore.point(:,end),xg_new);
        xstore.secant = v;
        relaxation_param_new = relaxation_param;
        xstore.point(:,end) = [];
    else
        slope = secant(xstore.point);
        xg = xstore.point(:,end) + ds*slope;
        [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param);
        v = secant(xstore.point(:,end),xg_new);
        while ~success && (ds>step_min)
            ds = ds/2;
            relaxation_param_new = relaxation_param_new*0.9;
            xg = xstore.point(:,end) + ds*slope;
            [xg_new,success,guess] = cont4(ds,xstore,Jx,f,xg,relaxation_param_new);
            v = secant(xstore.point(:,end),xg_new);
        end
        relaxation_param_new = relaxation_param;
    end
    V = xstore.secant(:,end);
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
    if ~success
        ds = -ds/2;
        xg = xstore.point(:,end) + ds*slope;
        [xg_new,success,guess] = cont4(-ds,xstore,Jx,f,xg,relaxation_param);
        check_v = secant(xstore.point);

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
% figure(1)
% for kk = 1:size(z.point,2)
%     hold on
%     text = sprintf("f_%i.jpg",kk);
%     scatter(z.point(3,kk),z.point(1,kk),'r')
%     axis([0, 40, 0, 5])
%     saveas(gcf,fullfile("pics",text))
%     % drawnow
%     % pause(1)
% end
% drawnow
% plot(z.point(3,:),z.point(1,:),'b')
% hold on
%
% figure(2)
% hold on
% plot(z.point(3,:),z.point(1,:),'b')
% for ii = 1:size(z.xg_sto,2)
% scatter(z.xg_sto{ii}(3,:),z.xg_sto{ii}(1,:))
% end
% scatter(z(1,end),z(2,end),'r','filled')
% x = linspace(min(z.point(1,:)),max(z.point(1,:)),100);
% plot(x,sin(x))



end