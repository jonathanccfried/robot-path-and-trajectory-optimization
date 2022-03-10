function [q1,q2,l] = gen3Rplanarpath(phi)

a1 =1;
a2 =1;
a3 =1;
l = 0:0.002:1;
l = l.^(2);
x = cos(2*pi*l.^2);
y = sin(2*pi*l.^2);

for i = 1:size(l,2)
   wx = x(i) - a3*cos(phi);
   wy = y(i) - a3*sin(phi);

   c2 = (wx^2 + wy^2 -a1^2 -a2^2)/(2*a1*a2);
   if c2 <=1
       s2_1 = sqrt(1 - c2^2);
       s2_2 = -sqrt(1-c2^2);
       q2_1 = atan2(s2_1,c2);
       q2_2 = atan2(s2_2,c2);

       s1_1 = (wy*(a1+a2*cos(q2_1))-a2*sin(q2_1)*wx)/(a1^2+a2^2+2*a1*a2*cos(q2_1));
       c1_1 = (wx*(a1+a2*cos(q2_1))+a2*sin(q2_1)*wy)/(a1^2+a2^2+2*a1*a2*cos(q2_1));

       s1_2 = (wy*(a1+a2*cos(q2_2))-a2*sin(q2_2)*wx)/(a1^2+a2^2+2*a1*a2*cos(q2_2));
       c1_2 = (wx*(a1+a2*cos(q2_2))+a2*sin(q2_2)*wy)/(a1^2+a2^2+2*a1*a2*cos(q2_2));
       q1_1 = atan2(s1_1,c1_1);
%        if q1_1 < 0
%            q1_1 = q1_1 + 2*pi;
%        end
        q1_2 = atan2(s1_2,c1_2);
%        if q1_2 < 0
%            q1_2 = q1_2 + 2*pi;
%        end

       q3_1 = phi - q1_1 -q2_1;
       q3_2 = phi - q1_2 -q2_2;

       q1(:,i) = [q1_1;q2_1;q3_1];
       q2(:,i) = [q1_2;q2_2;q3_2];
   end
   q1 = unwrap(q1,[],2);
   q2 = unwrap(q2,[],2);
end