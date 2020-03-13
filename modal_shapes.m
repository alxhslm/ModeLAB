function Vloc = modal_shapes(modes,geom)
Nmodes = size(modes.A,1);

%find least squared solution for each mode shape
Vloc = zeros(size(geom.r,1),Nmodes);
for j = 1:Nmodes    
    Aj = permute(modes.A(j,:,:),[3 2 1]);
    
    %find driving point hammer test with the largest response
    max_curr = -Inf;
    kMax = 1;
    for k = 1:length(geom.iAcc)
        iDrivePt = find(geom.bDrivePt(:,k) & geom.bModeAcc(k,j));
        if ~isempty(iDrivePt) && abs(Aj(k,iDrivePt)) > max_curr
             kMax = k;
             iMax = iDrivePt;
             max_curr = Aj(k,iDrivePt);
         end
    end
    
    fprintf('Using acc %d for mode %d\n',kMax,j)
    
    Phi = Aj ./ sqrt(Aj(kMax,iMax));
    Phi = abs(Phi) .*real(Phi);
 
    Vloc(:,j) = [Phi(kMax,geom.iFromHam) Phi(geom.iFromAcc,iMax).'];
end