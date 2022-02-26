function [maxlddoti, minlddoti,minjoint,maxjoint] = calc_inst_acc_cnstr(q,lambda,ldoti,qddotmax,qddotmin,step)

if step<3
qr = q(:,1:step+2);
lr = lambda(1:step+2);
gq = gradient(qr);
g2q = gradient(gq);
gl = gradient(lr);
g2l = gradient(gl);
dqdli = gq(:,step)./gl(step);
dqddli = g2q(:,step)./g2l(step);

elseif (3 <= step) && (step <= size(lambda,2)-2)
qr = q(:,step-2:step+2);
lr = lambda(step-2:step+2);
gq = gradient(qr);
g2q = gradient(gq);
gl = gradient(lr);
g2l = gradient(gl);
dqdli = gq(:,3)./gl(3);
dqddli = g2q(:,3)./g2l(3);

else
qr = q(:,step-2:end);
lr = lambda(step-2:end);
gq = gradient(qr);
g2q = gradient(gq);
gl = gradient(lr);
g2l = gradient(gl);
dqdli = gq(:,end-(size(lambda,2)-step))./gl(end-(size(lambda,2)-step));
dqddli = g2q(:,end-(size(lambda,2)-step))./g2l(end-(size(lambda,2)-step));
end


%Restrain Calc

lambdaqmax = (qddotmax - dqddli*ldoti^2)./dqdli;
lambdaqmin = (qddotmin - dqddli*ldoti^2)./dqdli;

maxlddoti = max([lambdaqmax,lambdaqmin],[],2);
minlddoti = min([lambdaqmax,lambdaqmin],[],2);

[maxlddoti,maxjoint] = min(maxlddoti);
[minlddoti,minjoint] = max(minlddoti);