function [theta,x,y] = eval_circle(p,x,y)
a = p(1);
b = p(2);
R = p(3);

theta = atan2(x-a,b-y);
y = b - R*cos(theta);
x = a + R*sin(theta);