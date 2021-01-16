function y = normalise(x,dim)
if nargin < 2
    dim = 1;
end
y = x ./L2norm(x,dim);