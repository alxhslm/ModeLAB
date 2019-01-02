function res = nyq_fit(w,H)
re = real(H);
im = imag(H);

bPlot = 0;
if bPlot
    fig = figure;
    axNyq = subplot(1,2,1);
    plot(re,im,'o')
    hold on
    axis equal
end

p = CircleFitByPratt([re im]);
R = p(3);
% p = CircleFitByAlex([re im]);
if bPlot
    plot_circle(p,'b');
end

theta = eval_circle(p,re,im);
theta = wrapToPi(theta);
wc = 0.5*(w(1:end-1) + w(2:end));
thetac = 0.5*(theta(1:end-1) + theta(2:end));

sweep = abs(diff(theta)./ (diff(w.^2)));

% coeff = polyfit(wc,sweep,3);
% wmax = roots(polyder(coeff));
% sweep_max = polyval(coeff,wmax);
% [sweep_max,iRes] = max(sweep);
% wmax = wc(iRes); 
% thetar = theta(iRes);

iFit = thetac>-2*pi/2 & thetac<2*pi/2;
coeff = polyfit(wc(iFit),thetac(iFit),3);
wmax = roots(polyder(polyder(coeff)));
sweep_max = interp1(wc(iFit),sweep(iFit),wmax,'linear','extrap');
% coeff = polyfit(wc(iFit),1./sweep(iFit),2);
% wmax = roots(polyder(coeff));
% sweep_max = 1./polyval(coeff,wmax);

% thetar = (round(mean(theta)/pi))*pi;
wr = wmax;%interp1(theta,w,thetar,'linear','extrap');
thetar = polyval(coeff,wmax);
%interp1(w,theta,wr,'linear','extrap');%polyval(polyint(coeff,theta(1)),wr);

% dtheta = theta - thetar;
% ii = find(dtheta(1:end-1).*dtheta(2:end) < 0,1);
% if ~isempty(ii)
%     wr = interp1(theta(ii:(ii+1)),w(ii:(ii+1)),thetar,'linear','extrap');
% else
%     wr = mean(w);
% end

if isnan(wr)
    keyboard
end

wh = interp1(theta,w,thetar+pi/2,'linear','extrap');
wl = interp1(theta,w,thetar-pi/2,'linear','extrap');
zr = (wh - wl)/2/wr;

if bPlot
    axAngle = subplot(2,2,2);
    plot(axAngle,w,theta)
    hold on
    plot(axAngle,wc,polyval(coeff,wc))
    plot(axAngle,wr,thetar,'o')

    axSweep = subplot(2,2,4);
    plot(axSweep,wc,1./sweep)
    hold on
    plot(axSweep,wmax,1./sweep_max,'o')
%     plot(axSweep,wc,polyval(coeff,wc),'r');
    
    plot(axNyq,p(1) + p(3)*sin(thetar),p(2) + p(3)*cos(thetar),'x')
%     plot(w,theta);
    waitforbuttonpress
    close(fig);
end

Ar = R * exp(1i*thetar);

res.wr = wr;
res.zr = zr;
res.Ar = Ar;

function plot_circle(p0,varargin)
theta = linspace(0,2*pi,100);
plot(p0(1)+p0(3)*sin(theta),p0(2)+p0(3)*cos(theta),varargin{:});

function t = myunwrap(t)
ii = 1;
while ~isempty(ii)
    ii = find(diff(t)<0,1);
    t((ii+1):end) = t((ii+1):end) + 2*pi;
end
% t = unwrap(t);