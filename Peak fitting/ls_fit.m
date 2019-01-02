function res = ls_fit(w,H)
bPlot = 0;

M = [H 2*1i*w.*H];
b = w.^2 .* H;
P = [real(M) 0*w-1 0*w;
     imag(M) 0*w 0*w-1];
q = [real(b);
     imag(b)];
c = P\q;

wr = sqrt(c(1));
zr = c(2)/wr;
Ar = c(3)+1i*c(4);

if bPlot
    H2 = Ar./(-w.^2 + 2*1i*w*wr*zr + wr^2);
    fig = figure;
    plot(w,abs(H),w,abs(H2))
    waitforbuttonpress
    close(fig)
end

res.wr = wr;
res.zr = zr;
res.Ar = Ar;