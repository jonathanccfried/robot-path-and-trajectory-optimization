function [ldotmin,t_final,qf,lf,indexmin] = constantpathspeedsolver(q,l,qdotmin,qdotmax,qddotmin,qddotmax,step)


    qf = [];
    lf = [];
for i =1:step:size(l,2)
    qf = [qf,q(:,i)];
    lf = [lf,l(i)];
end
dqdl = gradient(qf)./gradient(lf);
dqddl = [gradient(dqdl(1,:),lf);
        gradient(dqdl(2,:),lf);
        gradient(dqdl(3,:),lf);
        gradient(dqdl(4,:),lf);
        gradient(dqdl(5,:),lf);
        gradient(dqdl(6,:),lf)];

%% First Constraint Calculation

ldotmaxconstraint = zeros(size(dqdl));
for i = 1:size(dqdl,2)
ldotmaxconstraint(:,i) = max([qdotmax./dqdl(:,i),qdotmin./dqdl(:,i)],[],2);
end
sz1 = size(ldotmaxconstraint);
[ldotmaxconstraint,index1] = min(ldotmaxconstraint,[],'all');
[index1joint,index1path] = ind2sub(sz1,index1);
index1= [index1joint,index1path];
%% Second Constraint Calculation

ldotmaxconstraint2 = zeros(size(dqddl));
for i = 1:size(dqddl,2)
ldotmaxconstraint2(:,i) = sqrt(max([qddotmax./dqddl(:,i),qddotmin./dqddl(:,i)],[],2));
end
sz2 = size(ldotmaxconstraint2);
[ldotmaxconstraint2,index2] = min(ldotmaxconstraint2,[],'all');
[index2joint,index2path] = ind2sub(sz2,index2);
index2 = [index2joint,index2path];
%% 

ldotmin = min([ldotmaxconstraint2,ldotmaxconstraint]);
if ldotmaxconstraint < ldotmaxconstraint2
    indexmin = index1;
else
    indexmin = index2;
end
t_final = 1/ldotmin;

end