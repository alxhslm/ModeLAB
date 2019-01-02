function res = peak_fit(w,H)
N = length(H);
mag = abs(H);
ph = unwrap(angle(H));

%find the largest one
[magr,iRes] = max(mag);
if iRes > 1 && iRes < N && ((magr - mag(iRes-1)) * (mag(iRes+1) - magr)) < 0
    %then its a peak
    wr = w(iRes);
else
    iRes = ceil(N/2);
    wr = NaN;
end

iL = find(mag(1:iRes) < mag(iRes)/sqrt(2),1,'last');
iH = iRes + find(mag(iRes:end) < mag(iRes)/sqrt(2),1,'first')-1;

if isempty(iL) || isempty(iH)
    if isempty(iL)
        iL = 1;
    end
    if isempty(iH)
        iH = length(mag);
    end
    zr = (w(iH) - w(iL))/wr/(tan(abs(ph(iH)-ph(iRes))/2) + tan(abs(ph(iL)-ph(iRes))/2));
else
    zr = (w(iH) - w(iL))/wr/2;
end

Ar = H(iRes) * (2*1i*zr*wr^2);

res.wr = wr;
res.zr = zr;
res.Ar = Ar;