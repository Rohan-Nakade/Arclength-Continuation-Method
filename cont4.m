function [xg,successful,xg_store, A] = cont4(ds,xstore,Jx,f,xg,a)
    xg_store = xg;
    xn = xstore.point(:,end);
    dx = xg - xn;

    tol = norm(dx)*1e-2;
    maxiter = 2e3;
    numiter = 0;

%     theta = linspace(0,0.28,20);
%     theta2 = linspace(0,-pi/6,20);
%     figure(1)
%     hold on
%     plot(sin(theta),cos(theta),"blue","LineWidth",2)
%     plot(xn(2) + ds*cos(theta2),xn(1)+ds*sin(theta2),"green","LineWidth",2)

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

        Ds = norm(xg - xn);
        Dx = xg - xn;
        xg = xn + Dx*ds/Ds;
        
        numiter = numiter + 1;
        xg_store = [xg_store,xg];
%         figure(1)
%         hold on
%         if numiter == 1
%             line([xn(end),xg(end)],[xn(1),xg(1)],"Color","cyan","LineWidth",2)
%         end
%         scatter(xn(end),xn(1),100,"filled","diamond","MarkerFaceColor","#fb8500")
%         scatter(xg(end),xg(1),90,"filled")
%         ylim([0.96 1.015])
%         xlabel("x","FontSize",20,"Interpreter","latex")
%         ylabel("y","FontSize",20,"Interpreter","latex")
%         set(gca,"FontSize",20,"TickLabelInterpreter","latex")
%         set(gcf,"Position",[488 156.2000 903.4000 586.4000])
%         exportgraphics(gca,"test.gif","Append",true)
    end
    if norm(dx)<=tol
        successful = 1;
    else
        successful = 0;
    end
    figure(1)
    hold on
    scatter(xn(end),xn(1),"filled","blue")
%     scatter(xg_store(end,:),xg_store(1,:))
    ylim([-0.1 3])
    xlim([-5 5])
    xlabel("$\mu_{0}$","FontSize",20,"Interpreter","latex")
    ylabel("A","FontSize",20,"Interpreter","latex")
    set(gca,"FontSize",20,"TickLabelInterpreter","latex")
    set(gcf,"Position",[488 156.2000 903.4000 586.4000])
    exportgraphics(gca,"subcritical.gif","Append",true)
end