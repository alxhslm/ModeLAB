function Vloc = modal_shapes(modes,geom)
Nmodes = size(modes.A,1);

%find least squared solution for each mode shape
Vloc = zeros(size(geom.r,1),Nmodes);
for j = 1:Nmodes    
    V_exp = real(permute(modes.A(j,:,:),[3 2 1]));
    
    %find driving point hammer test with the largest response
    max_curr = -Inf;
    kMax = 1;
    for k = 1:length(geom.iAcc)
        iDrivePt = find(geom.bDrivePt(:,k));
        if ~isempty(iDrivePt) && V_exp(k,iDrivePt) > max_curr
             kMax = k;
             iMax = iDrivePt;
             max_curr = V_exp(k,iDrivePt);
         end
    end
    
    fprintf('Using acc %d for mode %d\n',kMax,j)
    
    %this amplitude is then used to scake the mode shape
    X0 = [V_exp(kMax,geom.iFromHam) V_exp(geom.iFromAcc,iMax)].'; 
    X0 = X0 / sqrt(V_exp(kMax,iMax));
    [g0,V_mod0] = objfun(X0,V_exp,geom.iHam,geom.iAcc);
    
    options = optimoptions('fminunc','OptimalityTolerance',1E-12,'StepTolerance',1E-12,'MaxFunctionEvaluations',10000,'Display','Off');
    
    [X,err,flag] = fminunc(@(X)objfun(X,V_exp,geom.iHam,geom.iAcc),X0,options);
    [g,V_mod] = objfun(X,V_exp,geom.iHam,geom.iAcc);
    Phi = X;

    Vloc(:,j) = Phi;
end

function [g,V_mod] = objfun(X,V_exp,iHam,iAcc)
Phi = X;
Scale = 1;
V_mod = Phi(iAcc)*(Phi(iHam).') * Scale^2;
g = sum(abs(V_exp(:) - V_mod(:)).^2);