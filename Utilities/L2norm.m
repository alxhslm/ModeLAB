function y = L2norm(x,dim)
if nargin < 2
    dim = 1;
end
y = sqrt(sum(abs(x).^2,dim));