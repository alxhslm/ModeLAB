function Par = CircleFitByAlex(XY)

x = XY(:,1); y = XY(:,2);
a = (max(x) + min(x))/2;
b = (max(y) + min(y))/2;
Rx = (max(x) - min(x))/2;
Ry = (max(y) - min(y))/2;
R = (Rx+Ry)/2;
p0 = [a b R]';

Par = modelfit(x,y,@eval_circle,p0,[],[],y);

function y = eval_circle(p,x,y)
a = p(1);
b = p(2);
R = p(3);

theta = atan2(y-b,x-a);
y = b + R*sin(theta);