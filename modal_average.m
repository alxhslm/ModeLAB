function modes = modal_average(modes,exp,setup)
%if we have used peak-fitting, we need to average over the different
%estimates of natural frequency etc.

[NHam,Nmodes,NAccel] = size(modes.peak);
for j = 1:Nmodes
    wj = zeros(NHam,NAccel);
    zj = zeros(NHam,NAccel);
    Aj = zeros(NHam,NAccel);
    for i = 1:NHam
        for k = 1:NAccel
            wj(i,k) = modes.peak(i,j,k).wr;
            zj(i,k) = modes.peak(i,j,k).zr;
            Aj(i,k) = modes.peak(i,j,k).Ar;
        end
    end
    
    iBad = isnan(wj) | isnan(zj) |  wj > setup.wBand(j,2) | wj < setup.wBand(j,1) | zj < 0;
    iBand = exp.w(2:end-1) > setup.wBand(j,1) & exp.w(2:end-1) < setup.wBand(j,2);
    
    wj(iBad) = NaN;
    zj(iBad) = NaN;
    
    [w_accept,iwReject] = deleteoutliers(wj,0.5);
    modes.omega(j) = mean(w_accept);
    iBad = iBad | iwReject;
    
    [z_accept,izReject] = deleteoutliers(zj,0.5);
    modes.zeta(j) = mean(z_accept);
    iBad = iBad | izReject;
    
    Aj = Aj./(wj.^2.*zj) .* (modes.omega(j)^2 * modes.zeta(j));
    
    H = 1./(modes.omega(j)^2 + 2*1i*modes.zeta(j)*modes.omega(j)*exp.w(iBand) - exp.w(iBand).^2);
    
    %recompute modal participation factors for "bad" dof
    A = reshape(H \ reshape(exp.H(iBand,:,:),sum(iBand),[]),[NHam NAccel]);
    Aj(iBad) = A(iBad);
    
    modes.A(j,:,:) = Aj;
end