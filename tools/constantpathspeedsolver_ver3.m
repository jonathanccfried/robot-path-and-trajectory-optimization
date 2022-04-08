function [ldotmin,t_final,indexmin,locktype,lag] = constantpathspeedsolver_ver3(qfit,qprime,qdoubleprime,l,qdotmin,qdotmax,qddotmin,qddotmax)


    qf = qfit;
    lf = l;
    qprimef = qprime;
    qdoubleprimef = qdoubleprime;

%% First Constraint Calculation

ldotmaxconstraint = zeros(size(qprimef));
for i = 1:size(qprimef,2)
ldotmaxconstraint(:,i) = max([qdotmax./qprimef(:,i),qdotmin./qprimef(:,i)],[],2);
end
sz1 = size(ldotmaxconstraint);
[ldotmaxconstraint,index1] = min(ldotmaxconstraint,[],'all');
[index1joint,index1path] = ind2sub(sz1,index1);
index1= [index1joint,index1path];
%% Second Constraint Calculation

ldotmaxconstraint2 = zeros(size(qdoubleprimef));
for i = 1:size(qdoubleprimef,2)
ldotmaxconstraint2(:,i) = sqrt(max([qddotmax./qdoubleprimef(:,i),qddotmin./qdoubleprimef(:,i)],[],2));
end
sz2 = size(ldotmaxconstraint2);
[ldotmaxconstraint2,index2] = min(ldotmaxconstraint2,[],'all');
[index2joint,index2path] = ind2sub(sz2,index2);
index2 = [index2joint,index2path];
%% 

ldotmin = min([ldotmaxconstraint2,ldotmaxconstraint]);
if ldotmaxconstraint < ldotmaxconstraint2
    indexmin = index1;
    locktype = ['vel'];
    ctrmax = qprimef(index1(1),index1(2))*ldotmin - qdotmax(index1(1));
    ctrmin = -qprimef(index1(1),index1(2))*ldotmin + qdotmin(index1(1));
    if ctrmin > ctrmax
        locktype = [locktype,'min'];
        lag = -(1/(ldotmin^2))*(1/qprimef(index1(1),index1(2)));
    else
        locktype = [locktype,'max'];
        lag = -(1/(ldotmin^2))*(-1/qprimef(index1(1),index1(2)));
    end    

else
    indexmin = index2;
    locktype = ['acc'];
    ctrmax = qdoubleprimef(index2(1),index2(2))*ldotmin^2 - qddotmax(index2(1));
    ctrmin = -qdoubleprimef(index2(1),index2(2))*ldotmin^2 + qddotmin(index2(1));
    if ctrmin > ctrmax
        locktype = [locktype,'min'];
        lag = -(1/(2*ldotmin^3))*(1/qdoubleprimef(index2(1),index2(2)));
    else
        locktype = [locktype,'max'];
        lag = -(1/(2*ldotmin^3))*(-1/qdoubleprimef(index2(1),index2(2)));
    end
end
t_final = 1/ldotmin;

end