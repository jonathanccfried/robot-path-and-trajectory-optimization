function [lambdaddothist,maxl2dothist,minl2dothist] = plot_l2dot_lddot(qf,lf,qddotmax,qddotmin,lddotmax,lddotmin,ldotmaxconstraint,ldotminconstraint,slice)

%slice gives a number of iterations k within the path for which a 2D plot
%is wanted

dqdl = gradient(qf)./gradient(lf);
dqddl = [gradient(dqdl(1,:),lf);
        gradient(dqdl(2,:),lf);
        gradient(dqdl(3,:),lf);
        gradient(dqdl(4,:),lf);
        gradient(dqdl(5,:),lf);
        gradient(dqdl(6,:),lf)];
samples = size(lf,2);


for n=1:samples
lambdaddotstep =  lddotmin(n):0.001:lddotmax(n);
minl2dotlist = zeros(1,size(lambdaddotstep,2));
maxl2dotlist = zeros(1,size(lambdaddotstep,2));
l2dotmax = ldotmaxconstraint(n)^2; %maximum path velocity defined by qdot constraints
l2dotmin = ldotminconstraint(n)^2; %minimum path velocity defined by qdot constraints

%Calculate ldot^2 for each lddot in lambdaddotstep
for k = 1:size(lambdaddotstep,2)
    dqdl(:,n);
    dqddl(:,n);
    l2dotqmax = (qddotmax - dqdl(:,n)*lambdaddotstep(k))./dqddl(:,n);
    l2dotqmin = (qddotmin - dqdl(:,n)*lambdaddotstep(k))./dqddl(:,n);

    l2dotmaxT = max([l2dotqmax,l2dotqmin],[],2); %velocity constraints that come from qddot constraints
    l2dotminT = min([l2dotqmax,l2dotqmin],[],2);

    [l2dotmaxT,indexmin] = min(l2dotmaxT);
    [l2dotminT,indexmax] = max(l2dotminT);

%     if n==83
%         if lambdaddotstep(k) <= -0.46303
%         indexmin
%         indexmax
%         end
%     end
    minl2dotlist(k) = max(l2dotminT,l2dotmin); %choose between velocity constraint given by qdot and qddot
    maxl2dotlist(k) = min(l2dotmaxT,l2dotmax);
end
if ismember(n,slice)
    figure(n)
    plot(lambdaddotstep,maxl2dotlist,'b',lambdaddotstep,minl2dotlist,'r')
    grid on
    title(['$\dot{\lambda}^2 \times \ddot{\lambda}$\,\mbox{for step }n=',num2str(n)],'Interpreter','latex')
xlabel('$\ddot{\lambda}$','Interpreter','latex')
ylabel('$\dot{\lambda}^2$','Interpreter','latex')
xlim([lddotmin(n)-1,lddotmax(n)+1]);
legend('$\dot{\lambda}^2_{max}\;\mbox{constr}$','$\dot{\lambda}^2_{min}\;\mbox{constr}$','Interpreter','latex')

end
lambdaddothist{n} = lambdaddotstep;
maxl2dothist{n} = maxl2dotlist;
minl2dothist{n} = minl2dotlist;
end
figure(n+1)
for n = 1:samples
    plot3(cell2mat(lambdaddothist(n)),cell2mat(maxl2dothist(n)),lf(n)*ones(size(cell2mat(lambdaddothist(n)))),'b',cell2mat(lambdaddothist(n)),cell2mat(minl2dothist(n)),lf(n)*ones(size(cell2mat(lambdaddothist(n)))),'r')
    hold on
end
grid on
    title(['$\lambda \times \dot{\lambda}^2 \times \ddot{\lambda}$ feasibility map'],'Interpreter','latex')
xlabel('$\ddot{\lambda}$','Interpreter','latex')
ylabel('$\dot{\lambda}^2$','Interpreter','latex')
zlabel('$\lambda$','Interpreter','latex')
hold off

figure(n+2)
for n = 1:samples
    plot(cell2mat(lambdaddothist(n)),lf(n)*ones(size(cell2mat(lambdaddothist(n)))),'b',cell2mat(lambdaddothist(n)),lf(n)*ones(size(cell2mat(lambdaddothist(n)))),'r')
    hold on
end
grid on
    title(['$\lambda \times \ddot{\lambda}$ feasibility map'],'Interpreter','latex')
xlabel('$\ddot{\lambda}$','Interpreter','latex')
ylabel('$\lambda$','Interpreter','latex')
hold off

figure(n+3)
for n = 1:samples
    plot(lf(n)*ones(size(cell2mat(lambdaddothist(n)))),cell2mat(maxl2dothist(n)),'b',lf(n)*ones(size(cell2mat(lambdaddothist(n)))),cell2mat(minl2dothist(n)),'r')
    hold on
end
grid on
    title(['$\lambda \times \dot{\lambda}^2$ feasibility map'],'Interpreter','latex')
ylabel('$\dot{\lambda}^2$','Interpreter','latex')
xlabel('$\lambda$','Interpreter','latex')
hold off
