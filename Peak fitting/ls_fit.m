function res = ls_fit(w,h)
bPlot = 0;

dims = size(h);
if length(dims)<3
    dims(end+1) = 1;
end
NDof = prod(dims(2:3));
H = h(:);
W = repmat(w,NDof,1);

M = [H 2*1i*W.*H];
b = W.^2 .* H;
P = [real(M) kron(eye(NDof),0*w-1) kron(zeros(NDof),0*w);
     imag(M) kron(zeros(NDof),0*w) kron(eye(NDof),0*w-1)];
q = [real(b);
     imag(b)];
c = P\q;

wr = sqrt(c(1));
zr = c(2)/wr;
Ar = c(2+(1:NDof))+1i*c(2+NDof+(1:NDof));

if bPlot
    h2 = Ar.'./(-w.^2 + 2*1i*w*wr*zr + wr^2);
    fig = figure;
    plot(w,abs(reshape(h,[],NDof)),w,abs(h2))
    waitforbuttonpress
    close(fig)
end

% figure
% plot(real(Ar),imag(Ar),'.')
% axis equal
% close all

res.wr = wr;
res.zr = zr;
res.Ar = reshape(Ar,dims(2:3));