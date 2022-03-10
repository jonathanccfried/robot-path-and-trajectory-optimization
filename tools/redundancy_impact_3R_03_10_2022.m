clear

phistep = -pi/2:0.01*pi:pi/2;
qdotmax = [1.75;1.57;1]; % joint speed limits in rad/s
qdotmin = -qdotmax;
qddotmax = 20*qdotmax; %joint acceleration limits in rad/s^2. Reasonable guess that it takes 0.1 seconds for a joint to reach its max limit
qddotmin = -qddotmax;

for i=1:size(phistep,2)
    [q1,q2,l] = gen3Rplanarpath(phistep(i));
    qlist(6*(i-1)+1:6*(i-1)+3,:) = q1;
    qlist(6*(i-1)+4:6*(i-1)+6,:) = q2;
    [ldotmin1,t_final1,~,~,indexmin1] = constantpathspeedsolver(q1,l,qdotmin,qdotmax,qddotmin,qddotmax,1);
    [ldotmin2,t_final2,~,~,indexmin2] = constantpathspeedsolver(q2,l,qdotmin,qdotmax,qddotmin,qddotmax,1);
    ldotminlist(2*i-1)=ldotmin1;
    ldotminlist(2*i) = ldotmin2;
    tlist(2*i-1) = t_final1;
    tlist(2*i) = t_final2;
    indexminlist(2*i-1,:) = indexmin1;
    indexminlist(2*i,:) = indexmin2;
end

[tmin,mindex]= min(tlist);
[tmax,maxdex] = max(tlist);
bestdex = indexminlist(mindex,:);
worstdex = indexminlist(maxdex,:);
bestlblock = l(bestdex(2));
worstlblock=l(worstdex(2));
bestjointblock=worstdex(1);
worstjointblock=worstdex(1);
ldotbest = ldotminlist(mindex);
ldotworst = ldotminlist(maxdex);
qbest = qlist(3*(mindex-1)+1:3*(mindex-1)+3,:);
qworst = qlist(3*(maxdex-1)+1:3*(maxdex-1)+3,:);
xbest = cos(qbest(1,:)) + cos(qbest(1,:)+qbest(2,:)) + cos(qbest(1,:) + qbest(2,:) + qbest(3,:));
ybest = sin(qbest(1,:)) + sin(qbest(1,:)+qbest(2,:)) + sin(qbest(1,:) + qbest(2,:) + qbest(3,:));
xworst = cos(qworst(1,:)) + cos(qworst(1,:)+qworst(2,:)) + cos(qworst(1,:) + qworst(2,:) + qworst(3,:));
yworst = sin(qworst(1,:)) + sin(qworst(1,:)+qworst(2,:)) + sin(qworst(1,:) + qworst(2,:) + qworst(3,:));
figure()
plot(l,qbest)
grid on
title(['$\lambda \times q$ with highest constant path speed'],'Interpreter','latex')
ylabel('$q$ (rad)','Interpreter','latex')
xlabel('$\lambda$','Interpreter','latex')
legend('$q_1$','$q_2$','$q_3$','Interpreter','latex')
figure()
plot(l,qworst)
grid on
title(['$\lambda \times q$ with lowest constant path speed'],'Interpreter','latex')
ylabel('$q$ (rad)','Interpreter','latex')
xlabel('$\lambda$','Interpreter','latex')
legend('$q_1$','$q_2$','$q_3$','Interpreter','latex')
figure()
plot(xworst,yworst)
grid on
title(['Cartesian Path with lowest constant path speed'],'Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
xlabel('$x$ m','Interpreter','latex')
figure()
plot(xbest,ybest)
grid on
title(['Cartesian Path with highest constant path speed'],'Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
xlabel('$x$ m','Interpreter','latex')
