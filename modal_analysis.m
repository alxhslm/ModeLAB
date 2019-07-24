function [modes,han] = modal_analysis(P,V,res_func)
% MODAL_ANALYSIS Performs modal analysis on multiple hammer tests
%
%      Inputs:
%           - P: Structure containing test configuration, such as locations 
%                of accelerometers, and natural frequencies bounds.
%           - V: Structure array containing response to each hammer test
%           - res_func: Method used to extract modal properties from each
%                peak. This can be 'peak_fit','nyq_fit','line_fit' or 'ls_fit'
%
%     Outputs:
%           - modes: Structure containig the modal properties and mode shapes       
%           - han: handles to the figures and axes used to plot the FRFs.
%                  Useful if you want to overlay information later.

if nargin < 3
    res_func = 'peak_fit';
end

if ~isfield(P,'fPhZero')
    P.fPhZero = 0;
end
% if ~isfield(P,'iSurface')
%     P.iSurface = 0*P.xHammer;
% end

NAccel = length(P.xAccel);
iAccel = zeros(NAccel,1);
nAcc = [P.axAccel;
        P.ayAccel;
        P.azAccel];
nHam = [P.FxHammer;
        P.FyHammer;
        P.FzHammer];
rHam = [P.xHammer;
        P.yHammer;
        P.zHammer];
rAcc = [P.xAccel;
        P.yAccel;
        P.zAccel];
   
        
%flip sign if accelerometer in the negative direction
sHam = (1-2*any(nHam < 0,1));
sAcc = (1-2*any(nAcc < 0,1));
sign = sHam'*sAcc;

nHam = nHam .* sHam;
nAcc = nAcc .* sAcc;

for i = 1:NAccel
    ii = find(abs(nAcc(:,i)'*nHam) > 0 & P.xHammer == P.xAccel(i)  & P.yHammer == P.yAccel(i) & P.zHammer == P.zAccel(i));
    if ~isempty(ii)
        iAccel(i) = ii;
    else
        iAccel(i) = NaN;
    end
end

%frequency bands containing each mode
wBand = [P.fL' P.fH']*2*pi;
Nmodes = size(wBand,1);

han.exp.fig = figure;
for k = 1:NAccel
    if P.axAccel(k) ~= 0
        label = 'Horizontal';
    elseif P.ayAccel(k) ~= 0
        label = 'Vertical';
    else
        label = 'Axial';
    end
        
    han.exp.axMag(k) = subplot(2,NAccel,0+k);
    hold on
    if k == 1
        ylabel(han.exp.axMag(k),'Mag')
    end
    set(han.exp.axMag(k),'yscale','log')
    title(label)
    
    han.exp.axPh(k) = subplot(2,NAccel,NAccel+k);  
    hold on
    if k == 1
        ylabel(han.exp.axPh(k),'Ph (deg)')
    end
    xlabel('f (Hz)')
end

Ntest = length(P.zHammer);
NSig = length(P.zAccel);

if size(V,1) > Ntest
    V = V(1:Ntest,:);
end

if size(V,1) > NSig
    V = V(:,1:NSig);
end

V = process_responses(V);

test_col = lines(Ntest);
mode_col = lines(Nmodes);

%convert acc/vel to disp
for i = 1:NSig
    w = V(1,i).Frequency + eps;
    if strncmp(V(1,i).Units,['m/s' char(178)],4)
        fprintf('Changing units of channel %d from ''m/s2'' to ''m''\n',i)
        scale = w.^2;
        units = ['m' V(1,i).Units(5:end)];
    elseif strncmp(V(1,i).Units,'m/s',3)
        fprintf('Changing units of channel %d from ''m/s'' to ''m''\n',i)
        scale = 1i*w;
        units = ['m' V(1,i).Units(4:end)];
    else
        scale = 0*w + 1;
        units = V(1,i).Units;
    end
            
    for j = 1:Ntest
        %flip sign if hammer in the negative direction

        %apply scaling
        V(j,i).H = sign(j,i)*V(j,i).H ./ scale;
        V(j,i).Real = real(V(j,i).H);
        V(j,i).Imaginary = imag(V(j,i).H);
        V(j,i).Magnitude = abs(V(j,i).H);
        V(j,i).Phase = phase(V(j,i).H);
        V(1,i).Units = units;
    end
end

disp('Plotting FRFs...')
for i = 1:Ntest
    for k = 1:NSig
        [~,ii] = min(abs(V(i,k).Frequency - P.fPhZero*2*pi));
        V(i,k).Phase = V(i,k).Phase - floor(V(i,k).Phase(ii));
        han.exp.hMag(i,k) = plot(han.exp.axMag(k),V(i,k).Frequency/2/pi,V(i,k).Magnitude,'-','color',test_col(i,:));
        han.exp.hPh(i,k)  = plot(han.exp.axPh(k) ,V(i,k).Frequency/2/pi,V(i,k).Phase*180/pi ,'color',test_col(i,:));
        
        if i == iAccel(k)
            set(han.exp.hMag(i,k),'LineWidth',2)
            set(han.exp.hPh(i,k),'LineWidth',2)
        end
    end
end

%highlight each mode
for j = 1:Nmodes
    wPlot = wBand(j,[1 2 2 1 1]);
    
    for k = 1:NAccel
        yPlot = ylim(han.exp.axMag(k)); yPlot = yPlot([1 1 2 2 1]);
        han.exp.hPatchMag(j,k) = patch(han.exp.axMag(k),wPlot/2/pi,yPlot,mode_col(j,:),'FaceAlpha',0.3,'EdgeColor','none');
        
        yPlot = ylim(han.exp.axPh(k)); yPlot = yPlot([1 1 2 2 1]);
        han.exp.hPatchPh(j,k)  = patch(han.exp.axPh(k) ,wPlot/2/pi,yPlot,mode_col(j,:),'FaceAlpha',0.3,'EdgeColor','none');
    end
end

%plot some horizontal lines every 180deg
for k = 1:NSig
    for i = 1:20
        han.exp.h180(i) = plot(han.exp.axPh(k),[0 max(V(1,k).Frequency)]/2/pi,-pi*i*[1 1]*180/pi,'k--');
    end
end

leg = {V(:,1).Label};
linkaxes([han.exp.axMag(:);han.exp.axPh(:)],'x');
drawnow

%% Now extract the frequencies and damping ratios
disp('Isolating peaks...')
if strcmp(res_func,'rfp')
    omega = V(1).Frequency;
    iFit = omega > 50 & omega < 2500;%100*2*pi & omega < 120*2*pi;
    nExtra = 1;
    nModes = Nmodes;
    wt = 0*omega-10;
    for j = 1:Nmodes
        %find the peaks in the specified range
        iBand = omega(2:end-1) > wBand(j,1) & omega(2:end-1) < wBand(j,2);
%         wt(iBand) = 1;
    end
    for i = 1:Ntest
        for k = 1:NSig
            rec = V(i,k).H;      
            modal_par = rfp(rec(iFit),omega(iFit),nModes+nExtra,wt(iFit));
            plot(han.exp.axMag(k),omega(iFit)/2/pi,abs(modal_par.alpha),'color',test_col(i,:),'LineStyle','--')
            plot(han.exp.axPh(k) ,omega(iFit)/2/pi,unwrap(angle(modal_par.alpha))*180/pi,'color',test_col(i,:),'LineStyle','--')
            for j = 1:nModes
                nyq(i,j,k).Frequency = modal_par.Frequency(j+nExtra);
                nyq(i,j,k).Damping = modal_par.Damping(j+nExtra);
                nyq(i,j,k).Magnitude = modal_par.Magnitude(j+nExtra);
                nyq(i,j,k).Phase = modal_par.Phase(j+nExtra);
                nyq(i,j,k).Real = modal_par.Real(j+nExtra);
                nyq(i,j,k).Imaginary = modal_par.Imaginary(j+nExtra);
            end
        end
    end
elseif strcmp(res_func,'grfp')
    omega = V(1).Frequency;
    iFit = omega > 50 & omega < 2500;%100*2*pi & omega < 120*2*pi;
    nExtra = 1;
    nModes = Nmodes;
    wt = 0*omega+1;
    for j = 1:nModes
        %find the peaks in the specified range
        iBand = omega(2:end-1) > wBand(j,1) & omega(2:end-1) < wBand(j,2);
        wt(iBand) = 1;
    end
    
    for k = 1:NSig
        for i = 1:Ntest
            rec(:,i+(k-1)*Ntest) = V(i,k).H;
%             modal_par = grfp(rec(iFit,i+(k-1)*Ntest),omega(iFit),nModes+nExtra,wt(iFit));
%             plot(han.exp.axMag(k),omega(iFit)/2/pi,abs(modal_par.alpha),'color',test_col(i,:),'LineStyle','--')
%             plot(han.exp.axPh(k) ,omega(iFit)/2/pi,unwrap(angle(modal_par.alpha))*180/pi,'color',test_col(i,:),'LineStyle','--')
%             for j = 1:nModes
%                 nyq(i,j,k).Frequency = modal_par.Frequency(j+nExtra);
%                 nyq(i,j,k).Damping = modal_par.Damping(j+nExtra);
%                 nyq(i,j,k).Magnitude = modal_par.Magnitude(j+nExtra);
%                 nyq(i,j,k).Phase = modal_par.Phase(j+nExtra);
%                 nyq(i,j,k).Real = modal_par.Real(j+nExtra);
%                 nyq(i,j,k).Imaginary = modal_par.Imaginary(j+nExtra);
%             end
        end
    end
    modal_par = grfp(rec(iFit,:),omega(iFit),nModes+nExtra,wt(iFit));

    for k = 1:NSig
        for i = 1:Ntest
            plot(han.exp.axMag(k),omega(iFit)/2/pi,abs(modal_par.alpha(:,i+(k-1)*Ntest)),'color',test_col(i,:),'LineStyle','--')
            plot(han.exp.axPh(k) ,omega(iFit)/2/pi,unwrap(angle(modal_par.alpha(:,i+(k-1)*Ntest)))*180/pi,'color',test_col(i,:),'LineStyle','--')
            for j = 1:nModes
                nyq(i,j,k).Frequency = modal_par.Frequency(j+nExtra);
                nyq(i,j,k).Damping = modal_par.Damping(j+nExtra);
                nyq(i,j,k).Magnitude = modal_par.Magnitude(i+(k-1)*Ntest,j+nExtra);
                nyq(i,j,k).Phase = modal_par.Phase(i+(k-1)*Ntest,j+nExtra);
                nyq(i,j,k).Real = modal_par.Real(i+(k-1)*Ntest,j+nExtra);
                nyq(i,j,k).Imaginary = modal_par.Imaginary(i+(k-1)*Ntest,j+nExtra);
            end
        end
    end
else
    for j = 1:Nmodes
        fprintf('- Mode %d\n',j)
        for i = 1:Ntest
            %node that the frequency base may be different between tests
            w = V(i,1).Frequency;

            %find the peaks in the specified range
            iBand = w(2:end-1) > wBand(j,1) & w(2:end-1) < wBand(j,2);

            for k = 1:NAccel
                %             dM_dw = adiff(V(i,k).Frequency,V(i,k).Magnitude);

                nyq(i,j,k) = feval(res_func,V(i,k).Frequency(iBand),V(i,k).H(iBand));
            end
        end
    end
end

disp('Extracting natural frequencies and modal damping...')
for i = 1:Ntest
    for j = 1:Nmodes
        for k = 1:NAccel
            W_exp(j,i,k) = nyq(i,j,k).wr;
            Z_exp(j,i,k) = nyq(i,j,k).zr;
            A_exp(j,i,k) = nyq(i,j,k).Ar;
        end
    end
end

for k = 1:NSig
    for i = 1:Ntest
        rec(:,i,k) = V(i,k).H;
    end
end

%     w_exp = mean(mean(W_exp,2,'omitnan'),3,'omitnan');
%     z_exp = mean(mean(Z_exp,2,'omitnan'),3,'omitnan');
w = V(1,1).Frequency;

for j = 1:Nmodes
    wj = squeeze(W_exp(j,:,:));
    zj = squeeze(Z_exp(j,:,:));
    Aj = squeeze(A_exp(j,:,:));
    
    iBad = isnan(wj) | isnan(zj) |  ...
            wj > wBand(j,2) | wj < wBand(j,1) | ...
            zj < 0;
    
    iBand = w(2:end-1) > wBand(j,1) & w(2:end-1) < wBand(j,2);
        
    wj(iBad) = NaN;
    zj(iBad) = NaN;
    
    [w_accept,iwReject] = deleteoutliers(wj,0.5);
    w_exp(j) = mean(w_accept);  
    iBad = iBad | iwReject;
    
    [z_accept,izReject] = deleteoutliers(zj,0.5);
    z_exp(j) = mean(z_accept);
    iBad = iBad | izReject;
    
    Aj = Aj./(wj.^2.*zj) .* (w_exp(j)^2 * z_exp(j));
    
    H = 1./(w_exp(j)^2 + 2*1i*z_exp(j)*w_exp(j)*w(iBand) - w(iBand).^2);
    
%     for k = 1:NSig
%         for i = 1:Ntest
%             b = [real(rec(iBand,i,k)); imag(rec(iBand,i,k))];
%             M = [real(H) -imag(H);
%                  imag(H)  real(H)];
%             x = M\b;
%             A1(i,k) = x(1) + 1i*x(2);
%         end
%     end

    %recompute modal participation factors for "bad" dof
    switch res_func
        case 'ls_fit'
            A = reshape(H \ reshape(rec(iBand,:,:),sum(iBand),[]),size(V));
        otherwise
            A = squeeze(interp1(w,rec,w_exp(j)))*(2*1i*w_exp(j)^2*z_exp(j));
    end

    Aj(iBad) = A(iBad);
    A_exp(j,:,:) = Aj;
end

%plot out the modal constants
for j = 1:Nmodes
    for i = 1:Ntest
        for k = 1:NAccel
            scale = (2*1i*w_exp(j)^2*z_exp(j));
            w = w_exp(j);
%             scale = (2*1i*W_exp(j,i,k)^2*Z_exp(j,i,k));
%             w = W_exp(j,i,k);
            Hmax = A_exp(j,i,k)./scale;
            plot(han.exp.axMag(k),w/2/pi,abs(Hmax),   'x','color',test_col(i,:));
            plot(han.exp.axPh(k) ,w/2/pi,angle(Hmax)*180/pi,'x','color',test_col(i,:));
        end
    end
end

%plot on FRF
for k = 1:NSig
    for j = 1:Nmodes
        han.exp.hResMag(j,k) = plot(han.exp.axMag(k),w_exp(j)/2/pi*[1 1],ylim(han.exp.axMag(k)),'color',mode_col(j,:));
        han.exp.hResPh(j,k)  = plot(han.exp.axPh(k) ,w_exp(j)/2/pi*[1 1],ylim(han.exp.axPh(k)) ,'color',mode_col(j,:));
    end
end
drawnow

%% Finally work out the mode shapes from the modal constants
disp('Extracting modeshapes...')
if ~isfield(P,'iSurface')
    iSurface = [];
else
    iSurface = P.iSurface;
end
[r_exp,n_exp,u_exp,iSurface,iHam,iAcc] = extract_modes(rAcc',nAcc',rHam',nHam',iSurface,A_exp);
[w_exp,iSort] = sort(w_exp);
z_exp  = z_exp(iSort);
u_exp  = u_exp(:,iSort);

modes.r = r_exp;
modes.n = n_exp;
modes.u = u_exp;
modes.zeta = z_exp;
modes.omega = w_exp;
modes.iSurf = iSurface;

w = V(1).Frequency';
U_exp = A_exp;
for k = 1:NSig
    for i = 1:Ntest
        U_exp(:,i,k) = u_exp(iHam(k),:).*u_exp(iAcc(k),:);
    end
end
% figure
for k = 1:NSig
    X = zeros(Ntest,length(w));
    for j = 1:Nmodes
        X = X + U_exp(j,:,k).'./ (-w.^2 + 2*1i*w_exp(j)*z_exp(j)*w + w_exp(j)^2);
    end
    for i = 1:Ntest
        u = X(i,:).';
        v = V(i,k).H; 
        
%         figure

%         axUser(1) = subplot(3,1,1);
%         plot(w/2/pi,imag(v),w/2/pi,imag(u))
%         hold on
%         semilogy(w/2/pi,abs(real(v)),w/2/pi,abs(real(u)))
%         hold on
%         
%         axUser(2) = subplot(3,1,2);
%         semilogy(w/2/pi,abs(imag(v)),w/2/pi,abs(imag(u)))
%         hold on
%         
%         axUser(3) = subplot(3,1,3);
%         semilogy(w/2/pi,abs(v),w/2/pi,abs(u))
%         hold on
        
        iFit = w > 50;
        if P.bResidMass && P.bResidStiffness
            [K(i,k),M(i,k)] = residual_mass_and_stiffness(w(iFit)',real(u(iFit)),real(v(iFit)));
        elseif P.bResidStiffness
            K(i,k) = residual_stiffness(w(iFit)',real(u(iFit)),real(v(iFit)));
            M(i,k) = Inf;
        elseif P.bResidMass
            K(i,k) = Inf;
            M(i,k) = residual_mass(w(iFit)',real(u(iFit)),real(v(iFit)));
        else
            K(i,k) = Inf;
            M(i,k) = Inf;
        end
            
        u = u - 1./(M(i,k)*w.^2)' + 1/K(i,k);
        
%         subplot(Ntest,NSig,(i-1)*NSig+k);
%         semilogy(w',abs([v u]))
           
        U(i,k).Frequency = w';
        U(i,k).H = u;
        U(i,k).Magnitude = abs(u);
        U(i,k).Phase = unwrap(angle(u));
        U(i,k).Real = real(u);
        U(i,k).Imaginary = imag(u);
        plot(han.exp.axMag(k),w/2/pi,U(i,k).Magnitude,'color',test_col(i,:),'LineStyle','--')
        plot(han.exp.axPh(k) ,w/2/pi,U(i,k).Phase*180/pi,'color',test_col(i,:),'LineStyle','--')
    end
end

disp('Plotting synthesised FRFs...')
han.model.fig = figure;
for k = 1:NAccel
    if P.axAccel(k) ~= 0
        label = 'Horizontal';
    elseif P.ayAccel(k) ~= 0
        label = 'Vertical';
    else
        label = 'Axial';
    end
        
    han.model.axMag(k) = subplot(2,NAccel,0+k);
    hold on
    if k == 1
        ylabel(han.model.axMag(k),'Mag')
    end
    set(han.model.axMag(k),'yscale','log')
    title(label)
    
    han.model.axPh(k) = subplot(2,NAccel,NAccel+k);  
    hold on
    if k == 1
        ylabel(han.model.axPh(k),'Ph (deg)')
    end
    xlabel('f (Hz)')
end

for i = 1:Ntest
    for k = 1:NSig
        han.model.hMag(i,k) = plot(han.model.axMag(k),U(i,k).Frequency/2/pi,U(i,k).Magnitude,'-','color',test_col(i,:));
        han.model.hPh(i,k)  = plot(han.model.axPh(k) ,U(i,k).Frequency/2/pi,U(i,k).Phase*180/pi ,'color',test_col(i,:));
        
        if i == iAccel(k)
            set(han.model.hMag(i,k),'LineWidth',2)
            set(han.model.hPh(i,k),'LineWidth',2)
        end
    end
end
linkaxes([han.model.axMag(:);han.model.axPh(:)],'x');

legend(han.exp.hMag(:,end),leg);
legend(han.model.hMag(:,end),leg);

Phi = u_exp;
modes.Phi = u_exp;
modes.M = Phi*Phi';
modes.C = Phi*diag(2*z_exp.*w_exp)*Phi';
modes.K = Phi*diag(w_exp)*Phi';
modes.Kres = K;
modes.Mres = M;

function [K,M] = residual_mass_and_stiffness(w,u,y)
c = [0*w+1 -1./(w.^2)]\(y-u);
K = 1/c(1);
M = 1/c(2);

function K = residual_stiffness(w,u,y)
K = 1/mean(y-u);

function M = residual_mass(w,u,y)
a = (-1./(w.^2))\(y-u);
M = -1./a;

function V = process_responses(V)
%checks whether we have time or frequency domain responses
%if in time domain, take FFT to yield the FRF

if ~isfield(V(1,1),'Freq') && isfield(V(1,1),'Time')
    %time domain
    NSig = size(V,2)-1;
    for i = 1:Ntest

        x = V(i,1).Real;
        t = V(i,1).Time;

        y = repmat(0*x,1,NSig);
        for k = 1:NSig
            y(:,k) = V(i,k+1).Real;
        end
        
        %remove nans
        iKeep = ~isnan(x) & ~any(isnan(y),2);
        x = x(iKeep);
        y = y(iKeep,:);
        t = t(iKeep);
        
%         ws = 2*pi./mean(diff(t));
%         wc = 200*2*pi;
%         plot(ax(1),t,y(:,1))
%         plot(ax(2),t,y(:,2))
        %apply butterworth filter
%         [num,den] = butter(2,wc/(ws/2));
%         y = filtfilt(num,den,y);
%         figure
%         plot(t,y)
%         for k = 1:size(y,2)
%             y(:,k) = smooth(y(:,k),20);
%         end
%         hold on
%         plot(t,y)
        
        %take fft
        X = fft(x);
        Y = fft(y);
        H = Y./X;

        ws = 2*pi./mean(diff(t));
        Nfft = length(t);
        w = ((1:Nfft)-1)'/Nfft * ws;
        
        for k = 1:NSig
            V2(i,k).Magnitude = abs(H(:,k));
            V2(i,k).Phase = angle(H(:,k));

            V2(i,k).Imaginary = imag(H(:,k));
            V2(i,k).Real = real(H(:,k));

            V2(i,k).Frequency = w;
        end
        
    end
    V = V2;
end