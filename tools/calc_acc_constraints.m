function [lambdadotminconst,lambdadotmaxconst,lambdaddotmin, lambdaddotmax,mindex,maxindex] = calc_acc_constraints(q,lambda,qdotmax,qdotmin,qddotmax,qddotmin,step)

% Objective > Find constraints lambdaddotmax and lambdaddotmin, given
% qddot = dqddlambda*dlambda2 + dqdlambda*lambdaddot
% step = plot lambdaddot x lambdadot on that particular step of path lambda

dqdlambda = gradient(q)./gradient(lambda);
dqddlambda = gradient(gradient(q))./gradient(gradient(lambda));
%dqdlambda2 = diff(q,1,2)./diff(lambda,1);
%dqddlambda2 = diff(q,2,2)./diff(lambda,2);

% First constraint qdotmin <= c'(lambda)*lambdadot <= qdotmax &&
% lambdadot > 0

lambdadotmaxconst = zeros(size(dqdlambda));
for i = 1:size(dqdlambda,2)
lambdadotmaxconst(:,i) = max([qdotmax./dqdlambda(:,i),qdotmin./dqdlambda(:,i)],[],2);
end
lambdadotmaxconst = min(lambdadotmaxconst,[],1);


lambdadotminconst = zeros(size(dqdlambda));
for i = 1:size(dqdlambda,2)
lambdadotminconst(:,i) = min([qdotmax./dqdlambda(:,i),qdotmin./dqdlambda(:,i)],[],2);
lambdadotminconst(:,i) = max([lambdadotminconst(:,i),zeros(size(lambdadotminconst(:,i)))],[],2);
end
lambdadotminconst = max(lambdadotminconst,[],1);


% Calculate range of lambdaddot

lambdaddotmin = zeros(size(q));
lambdaddotmax = zeros(size(q));

for i = 1:size(lambda,2)

% Calculate for each joint
size(q,1)
for j = 1:size(q,1)

%Check sign of dqdlambda and dqddlambda

if dqdlambda(j,i) == 0
    lambdaddotmax(j,i) = 0;
    lambdaddotmin(j,i) = 0;
elseif dqddlambda(j,i) >= 0 && dqdlambda(j,i) > 0

    lambdaddotmax(j,i) = (qddotmax(j)- dqddlambda(j,i)*lambdadotminconst(i)^2)/dqdlambda(j,i);
    lambdaddotmin(j,i) = (qddotmin(j) - dqddlambda(j,i)*lambdadotmaxconst(i)^2)/dqdlambda(j,i);
elseif dqddlambda(j,i) <= 0 && dqdlambda(j,i) < 0

    lambdaddotmax(j,i) = (qddotmin(j)- dqddlambda(j,i)*lambdadotminconst(i)^2)/dqdlambda(j,i);
    lambdaddotmin(j,i) = (qddotmax(j) - dqddlambda(j,i)*lambdadotmaxconst(i)^2)/dqdlambda(j,i);

elseif dqddlambda(j,i) <= 0 && dqdlambda(j,i) > 0

    lambdaddotmax(j,i) = (qddotmax(j)- dqddlambda(j,i)*lambdadotmaxconst(i)^2)/dqdlambda(j,i);
    lambdaddotmin(j,i) = (qddotmin(j) - dqddlambda(j,i)*lambdadotminconst(i)^2)/dqdlambda(j,i);
elseif dqddlambda(j,i) >= 0 && dqdlambda(j,i) < 0

    lambdaddotmax(j,i) = (qddotmin(j)- dqddlambda(j,i)*lambdadotmaxconst(i)^2)/dqdlambda(j,i);
    lambdaddotmin(j,i) = (qddotmax(j) - dqddlambda(j,i)*lambdadotminconst(i)^2)/dqdlambda(j,i);

end
end
end

%Choose minimum lambdaddotmax and maximum lambdaddotmin

[lambdaddotmax,mindex] = min(lambdaddotmax,[],1);
[lambdaddotmin,maxindex] = max(lambdaddotmin,[],1);

%Plot on step

n = step;
lambdaddotstep =  lambdaddotmin(n):0.001:lambdaddotmax(n);

%mindex 

if dqddlambda(mindex(n),n) >= 0 && dqdlambda(mindex(n),n) > 0
    ldot2max = (qddotmax(mindex(n)) -dqdlambda(mindex(n),n).*lambdaddotstep)./dqddlambda(mindex(n),n);
elseif dqddlambda(mindex(n),n) <= 0 && dqdlambda(mindex(n),n) < 0
    ldot2max = (qddotmin(mindex(n)) -dqdlambda(mindex(n),n).*lambdaddotstep)./dqddlambda(mindex(n),n);
elseif dqddlambda(mindex(n),n) <= 0 && dqdlambda(mindex(n),n) > 0
    ldot2max = (qddotmax(mindex(n)) -dqdlambda(mindex(n),n).*lambdaddotstep)./dqddlambda(mindex(n),n);   
elseif dqddlambda(mindex(n),n) >= 0 && dqdlambda(mindex(n),n) < 0
    ldot2max = (qddotmin(mindex(n)) -dqdlambda(mindex(n),n).*lambdaddotstep)./dqddlambda(mindex(n),n);

end

%maxindex calc

if dqddlambda(maxindex(n),n) >= 0 && dqdlambda(maxindex(n),n) > 0
    ldot2min = (qddotmin(maxindex(n)) -dqdlambda(maxindex(n),n).*lambdaddotstep)./dqddlambda(maxindex(n),n);
elseif dqddlambda(maxindex(n),n) <= 0 && dqdlambda(maxindex(n),n) < 0
    ldot2min = (qddotmax(maxindex(n)) -dqdlambda(maxindex(n),n).*lambdaddotstep)./dqddlambda(maxindex(n),n);
elseif dqddlambda(maxindex(n),n) <= 0 && dqdlambda(maxindex(n),n) > 0
    ldot2min = (qddotmin(maxindex(n)) -dqdlambda(maxindex(n),n).*lambdaddotstep)./dqddlambda(maxindex(n),n);
elseif dqddlambda(maxindex(n),n) >= 0 && dqdlambda(maxindex(n),n) < 0
    ldot2min = (qddotmax(maxindex(n)) -dqdlambda(maxindex(n),n).*lambdaddotstep)./dqddlambda(maxindex(n),n);
end

% Plots


%ldot2max = (qddotmax(mindex(n)) -dqdlambda(mindex(n),n).*lambdaddotstep)./dqddlambda(mindex(n),n);
%ldot2min = (qddotmin(maxindex(n)) -dqdlambda(maxindex(n),n).*lambdaddotstep)./dqddlambda(maxindex(n),n);
l2dotmin = lambdadotminconst(n)^2*ones(size(lambdaddotstep));
l2dotmax = lambdadotmaxconst(n)^2*ones(size(lambdaddotstep));

figure(1)
plot(lambdaddotstep,ldot2max,lambdaddotstep,ldot2min,lambdaddotstep,l2dotmin)
title('$\dot{\lambda}^2 \times \ddot{\lambda}$','Interpreter','latex')
xlabel('$\ddot{\lambda}$','Interpreter','latex')
ylabel('$\dot{\lambda}^2$','Interpreter','latex')
xlim([lambdaddotmin(n)-1,lambdaddotmax(n)+1])
legend('$\ddot{\lambda}_{max}\;\mbox{constr}$','$\ddot{\lambda}_{min}\;\mbox{constr}$','$\dot{\lambda} \geq 0$','Interpreter','latex')
figure(2)
plot(lambdaddotstep,ldot2max)
figure(3)
plot(lambdaddotstep,ldot2min)
figure(4)
plot(lambdaddotstep,l2dotmax)
figure(5)
plot(lambdaddotstep,l2dotmin)
end


