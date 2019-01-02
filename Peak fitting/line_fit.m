function peak = line_fit(w,H)
re = real(H);
im = imag(H);

mag = abs(H);
ph = angle(H);

iPivot = ceil(linspace(1,length(w),100));
O = w(iPivot);
for i = 1:length(iPivot)
    [mR(i,1),mI(i,1)] = line_fit_at(w,re,im,iPivot(i));
end

preal = polyfit(O.^2,mR,1);
pimag = polyfit(O.^2,mI,1);

% figure
% subplot(2,1,1)
% plot(O.^2,mR);
% hold on
% plot(O.^2,polyval(preal,O.^2))
% 
% subplot(2,1,2)
% plot(O.^2,mI);
% hold on
% plot(O.^2,polyval(pimag,O.^2))

nR = preal(1);
dR = preal(2);
nI = pimag(1);
dI = pimag(2);
p = nI/nR;
q = dI/dR;

eta = (q-p)/(1+p*q);
omega = sqrt(dR/((p*eta-1)*nR));
ar = omega^2*(p*eta-1)/((1+p^2)*dR);
br = -ar*p;
zeta = eta / 2;

Ar = ar + 1i*br;

peak.wr = omega;
peak.zr = zeta;
peak.Ar = Ar;

function [mR,mI] = line_fit_at(w,re,im,iCentre)
O = w(iCentre);

re_ = (re - re(iCentre))+eps;
im_ = (im - im(iCentre))+eps;

D = (w.^2 - O.^2)./(re_ + 1i*im_);

preal = polyfit(w.^2,real(D),1);
pimag = polyfit(w.^2,imag(D),1);

mR = preal(1);
mI = pimag(1);
cR = preal(2);
cI = pimag(2);

% figure
% subplot(1,2,1)
% plot(w.^2,real(D));
% hold on
% plot(w.^2,polyval(preal,w.^2));
% 
% subplot(1,2,2)
% plot(w.^2,imag(D));
% hold on
% plot(w.^2,polyval(pimag,w.^2));

