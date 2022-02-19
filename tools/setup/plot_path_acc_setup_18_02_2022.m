q = readmatrix('qcurve.txt'); %curve in joint space
q = q';
qdotmax = [1.75;1.57;1.57;2.97;2.09;3.32]; % joint speed limits in rad/s
qdotmin = -qdotmax;
qddotmax = [5;5;5;7;8;9]; %joint acceleration limits in rad/s^2
qddotmin = -qddotmax;
p = readmatrix('pcurve.txt'); % Curve in the cartesian frame
pcartesian = (p(:,1:3))'; % XYZ curve