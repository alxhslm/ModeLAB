function [rLoc,nLoc,Vloc,iSurface,iHam,iAcc] = extract_modes(rAcc,nAcc,rHam,nHam,iSurface,C_exp)

%identify how many unique points we have
rLoc = [rHam;
        rAcc];
    
nLoc = [nHam;
        nAcc]; 

A = [rLoc nLoc];
[C,iA,iC] = unique(A,'rows','stable');
rLoc = C(:,1:3);
nLoc = C(:,4:6);

iFromHam = iA(iA <= size(rHam,1));
iFromAcc = iA(iA > size(rHam,1)) - size(rHam,1);

NHam = size(rHam,1);
iHam = iC(1:NHam);
iAcc = iC(NHam+1:end);

if ~isempty(iSurface)
    iSurface = iSurface(:,iA)';
end

%find driving point responses
for k = 1:length(iAcc)
    ii = find(iHam == iAcc(k));
    if ~isempty(ii)
        iDrivF(k) = ii;
    else
        iDrivF(k) = NaN;
    end
end
Nmodes = size(C_exp,1);

%find least squared solution for each mode shape
for j = 1:Nmodes    
    V_exp = real(permute(C_exp(j,:,:),[3 2 1]));
    
    %find driving point hammer test with the largest response
    max_curr = -Inf;
    kMax = 1;
    for k = 1:length(iAcc)
         if ~isnan(iDrivF(k)) && V_exp(k,iDrivF(k)) > max_curr
             kMax = k;
             max_curr = V_exp(k,iDrivF(k));
         end
    end
    
    fprintf('Using acc %d for mode %d\n',kMax,j) 
    
    %this amplitude is then used to scake the mode shape
    X0 = [V_exp(kMax,iFromHam) V_exp(iFromAcc,iDrivF(kMax))].'; 
    X0 = X0 / sqrt(V_exp(kMax,iDrivF(kMax)));
%     Phi = abs(X0).*sign(real(X0));
%     Phi = imag(X0);
%     Phi = real(X0);
%     Scale = Phi\X0(iA);
%     X0 = [Phi; abs(Scale); angle(Scale)];
%     X0 = Phi;
    
    [g0,V_mod0] = objfun(X0,V_exp,iHam,iAcc);
    
    options = optimoptions('fminunc','OptimalityTolerance',1E-12,'StepTolerance',1E-12,'MaxFunctionEvaluations',10000,'Display','Off');
    
    [X,err,flag] = fminunc(@(X)objfun(X,V_exp,iHam,iAcc),X0,options);
    [g,V_mod] = objfun(X,V_exp,iHam,iAcc);
    Phi = X;
%     [V_exp.' V_mod.']
%     Phi = X(1:end-2);
    Vloc(:,j) = Phi;
end

function [g,V_mod] = objfun(X,V_exp,iHam,iAcc)
Phi = X;
Scale = 1;
% Phi = X(1:end-2);
% Scale = X(end-1)*exp(1i*X(end));

V_mod = Phi(iAcc)*(Phi(iHam).') * Scale^2;
g = sum(abs(V_exp(:) - V_mod(:)).^2);