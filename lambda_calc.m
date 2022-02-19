function lambda = lambda_calc(pcurve)

lambda = zeros(1,size(pcurve,2));
for i=1:(size(pcurve,2)-1)
    dist = pcurve(:,i+1) - pcurve(:,i);
    lambda(i+1) = norm(dist) + lambda(i);
end

lambda = lambda/(lambda(end));