clear

q = readmatrix('qcurve2.txt');
q = q';

p = readmatrix('pcurve2.txt');
pcartesian = (p(:,1:3))';
l = lambda_calc(pcartesian);
n = 15; %order of the poly

for i=1:size(q,1)
    a(i,:) = polyfit(l,q(i,:),n);
    qfit(i,:) = polyval(a(i,:),l);
    err(i) = norm(q(i,:)-qfit(i,:));
%     figure(i)
%     plot(l,q(i,:),l,polyval(a(i,:),l))
%     legend(strcat('q' , num2str(i)),'qfit')
    aprime(i,:) = polyder(a(i,:));
    adoubleprime(i,:) = polyder(aprime(i,:));
    qprime(i,:) = polyval(aprime(i,:),l);
    qdoubleprime(i,:) = polyval(adoubleprime(i,:),l);
end
err;

qdotmax = [1.75;1.57;1.57;2.97;2.09;3.32]; % joint speed limits in rad/s
qdotmin = -qdotmax;
qddotmax = (pi/180)*[312;292;418;2407;1547;3400]; %joint acceleration limits in rad/s^2. Reasonable guess that it takes 0.1 seconds for a joint to reach its max limit
qddotmin = -qddotmax;

tic
%[ldotmax, lag] = cvx_minimizer_lambda_constant_ver2(qfit,qprime,qdoubleprime,l,1,qdotmin,qdotmax,qddotmin,qddotmax);
[ldotmin,t_final,indexmin,locktype,lag] = constantpathspeedsolver_ver3(qfit,qprime,qdoubleprime,l,qdotmin,qdotmax,qddotmin,qddotmax);
toc