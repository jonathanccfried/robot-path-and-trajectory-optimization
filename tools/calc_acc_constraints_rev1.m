function [lf,qf,ldotminconstraint,ldotmaxconstraint,lddotmin, lddotmax,mindex,maxindex] = calc_acc_constraints_rev1(q,lambda,qdotmax,qdotmin,qddotmax,qddotmin,step)

% Objective > Find constraints lambdaddotmax and lambdaddotmin, given
% qddot = dqddlambda*dlambda2 + dqdlambda*lambdaddot
    qf = [];
    lf = [];
for i =1:step:size(lambda,2)
    qf = [qf,q(:,i)];
    lf = [lf,lambda(i)];
end
dqdl = gradient(qf)./gradient(lf);
dqddl = gradient(gradient(qf))./gradient(gradient(lf));
%dqdlambda2 = diff(q,1,2)./diff(lambda,1);
%dqddlambda2 = diff(q,2,2)./diff(lambda,2);

% First constraint qdotmin <= c'(lambda)*lambdadot <= qdotmax &&
% lambdadot > 0

ldotmaxconstraint = zeros(size(dqdl));
for i = 1:size(dqdl,2)
ldotmaxconstraint(:,i) = max([qdotmax./dqdl(:,i),qdotmin./dqdl(:,i)],[],2);
end
ldotmaxconstraint = min(ldotmaxconstraint,[],1);


ldotminconstraint = zeros(size(dqdl));
for i = 1:size(dqdl,2)
ldotminconstraint(:,i) = min([qdotmax./dqdl(:,i),qdotmin./dqdl(:,i)],[],2);
ldotminconstraint(:,i) = max([ldotminconstraint(:,i),zeros(size(ldotminconstraint(:,i)))],[],2);
end
ldotminconstraint = max(ldotminconstraint,[],1);


% Calculate range of lambdaddot

lddotmin = zeros(size(qf));
lddotmax = zeros(size(qf));

for i = 1:size(lf,2)
% Calculate for each joint


lmax1 = (qddotmax- dqddl(:,i).*ldotminconstraint(i)^2)./dqdl(:,i);
lmax2 = (qddotmax- dqddl(:,i).*ldotmaxconstraint(i)^2)./dqdl(:,i);
lmin1 = (qddotmin- dqddl(:,i).*ldotmaxconstraint(i)^2)./dqdl(:,i);
lmin2 = (qddotmin- dqddl(:,i).*ldotminconstraint(i)^2)./dqdl(:,i);



lddotmax(:,i) = max([lmax1,lmax2,lmin1,lmin2],[],2);
lddotmin(:,i) = min([lmax1,lmax2,lmin1,lmin2],[],2);
%lddotmax(:,i)
%lddotmin(:,i)
end

%Choose minimum lambdaddotmax and maximum lambdaddotmin among joints

[lddotmax,mindex] = min(lddotmax,[],1);
[lddotmin,maxindex] = max(lddotmin,[],1);

% Plot and check for feasibility


for n = 1:size(lf,2)

    lambdaddotstep =  lddotmin(n):0.01:lddotmax(n);
minl2dotlist = zeros(1,size(lambdaddotstep,2));
maxl2dotlist = zeros(1,size(lambdaddotstep,2));
l2dotmax = ldotmaxconstraint(n)^2; %maximum path velocity defined by qdot constraints
l2dotmin = ldotminconstraint(n)^2; %minimum path velocity defined by qdot constraints

%Calculate ldot^2 for each lddot in lambdaddotstep

for k = 1:size(lambdaddotstep,2)
    l2dotqmax = (qddotmax - dqdl(:,n)*lambdaddotstep(k))./dqddl(:,n);
    l2dotqmin = (qddotmin - dqdl(:,n)*lambdaddotstep(k))./dqddl(:,n);

    l2dotmaxT = max([l2dotqmax,l2dotqmin],[],2); %velocity constraints that come from qddot constraints
    l2dotminT = min([l2dotqmax,l2dotqmin],[],2);

    l2dotmaxT = min(l2dotmaxT);
    l2dotminT = max(l2dotminT);

    minl2dotlist(k) = max(l2dotminT,l2dotmin); %choose between velocity constraint given by qdot and qddot
    maxl2dotlist(k) = min(l2dotmaxT,l2dotmax);
end
    l2dotrange = abs(maxl2dotlist - minl2dotlist);
    %[~,zeroindex] = mink(l2dotrange,2);
%     if n == 240
%         lf(n)
%         plot(lambdaddotstep,maxl2dotlist,lambdaddotstep,minl2dotlist)
%     end

    %lddotmin(n) = min(lambdaddotstep(zeroindex(2)),lambdaddotstep(zeroindex(1)));
    %lddotmax(n) = max(lambdaddotstep(zeroindex(2)),lambdaddotstep(zeroindex(1)));
    BW = imregionalmin(l2dotrange);
    indexmin = find(BW);
    lddotmin(n) = lambdaddotstep(indexmin(1));
    lddotmax(n) = lambdaddotstep(indexmin(2));
end
end