function [modes,model,han] = modal_analysis(dataroot,options)
% MODAL_ANALYSIS Performs modal analysis on multiple hammer tests
%
%      Inputs:
%           - dataroot: Path to folder containing:
%                   * Setup CSV containing locations of hammers, accelerometers,
%                   natural frequencies bounds etc
%                   * ASCII FRF files from Signal Calc
%           - options: Structure constaining: 
%                   * 'method' : Method used to extract modal properties 
%                   * 'bResidMass' and 'bResidStiffness' : Residual mass/stiffness
%     Outputs:
%           - modes: Structure containing the modal properties and mode shapes       
%           - model: Structure containing the synthesised FRFs
%           - han: handles to the figures and axes used to plot the FRFs.
%                  Useful if you want to overlay information later.

if nargin < 2
    options = struct();
end

setup_csv_file = fullfile(dataroot, 'setup.csv'); 
setup_mat_file = fullfile(dataroot, 'setup.mat'); 
if ~isfile(setup_mat_file) || moddate(setup_csv_file) > moddate(setup_mat_file)
    setup = modal_setup(setup_csv_file);
    save(setup_mat_file,'-struct','setup')
else
    setup = load(setup_mat_file);
end

%% Extract FRFs and store in matrix
exp_mat_file = fullfile(dataroot, 'exp.mat'); 
exp = modal_load_frfs(dataroot);   

%% Trim frf / geomtery

if size(exp.H,2) > size(setup.rHam,1)
    %trim frf
    warning('Not enough hammer points specified in setup csv')
    exp.H = exp.H(:,1:size(setup.rHam,1),:);
    exp.TestLabel = exp.TestLabel(1:size(setup.rHam,1));
elseif size(exp.H,2) < size(setup.rHam,1)
    %trim setup
    warning('Too many hammer points specified in setup csv')
    setup.rHam = setup.rHam(1:size(exp.H,2),:);
    setup.nHam = setup.nHam(1:size(exp.H,2),:);
    setup.iBodyHam = setup.iBodyHam(1:size(exp.H,2));
    setup.sTest = setup.sTest(1:size(exp.H,2),:);
    setup.bDrivePt = setup.bDrivePt(1:size(exp.H,2),:);
end

if size(exp.H,3) > size(setup.rAcc,1)
    %trim frf
    warning('Not enough accelerometer points specified in setup csv')
    exp.H = exp.H(:,:,1:size(setup.rAcc,1));
elseif size(exp.H,3) < size(setup.rAcc,1)
    %trim setup
    warning('Too many accelerometer points specified in setup csv')
    setup.rAcc = setup.rAcc(1:size(exp.H,3),:);
    setup.nAcc = setup.nAcc(1:size(exp.H,3),:);
    setup.iBodyAcc = setup.iBodyAcc(1:size(exp.H,3));
    
    setup.sTest = setup.sTest(1:size(exp.H,3),:);
    setup.bDrivePt = setup.bDrivePt(1:size(exp.H,3),:);
    setup.AccLabel = setup.AccLabel(1:size(exp.H,3));
end

NAccel = size(exp.H,3);
NHam  = size(exp.H,2);
Nfreq = length(exp.w);

for k = 1:NAccel
    for i = 1:NHam
        %flip sign if hammer in the negative direction
        exp.H(:,i,k) = setup.sTest(i,k)*exp.H(:,i,k);
    end
end

Nmodes = size(setup.wBand,1);
options = default_options(options,NHam,Nmodes);

%% Plot FRFs
disp('Plotting FRFs...')
han.exp = modal_plot_frfs(exp,setup,options);

%highlight each mode
for j = 1:Nmodes
    wPlot = setup.wBand(j,[1 2 2 1 1]);
    
    for k = 1:NAccel
        yPlot = ylim(han.exp.axMag(k)); yPlot = yPlot([1 1 2 2 1]);
        han.exp.hPatchMag(j,k) = patch(han.exp.axMag(k),wPlot/2/pi,yPlot,options.mode_col(j,:),'FaceAlpha',0.3,'EdgeColor','none');
        
        yPlot = ylim(han.exp.axPh(k)); yPlot = yPlot([1 1 2 2 1]);
        han.exp.hPatchPh(j,k)  = patch(han.exp.axPh(k) ,wPlot/2/pi,yPlot,options.mode_col(j,:),'FaceAlpha',0.3,'EdgeColor','none');
    end
end

%% Extract natural frequencies and damping ratios
modes_mat_file = fullfile(dataroot, 'modes.mat'); 
if ~isfile(modes_mat_file) || moddate(setup_mat_file) > moddate(modes_mat_file) || moddate(exp_mat_file) > moddate(modes_mat_file)
    disp('Isolating peaks..')
    modes.fit.w = exp.w;
    modes.fit.H = NaN(Nfreq,NHam,NAccel);
    switch options.method
        case 'grfp'
            iFit = exp.w > setup.wMin & exp.w < setup.wMax;
            frf_mat = reshape(exp.H,length(exp.w),[]);
            modal_par = grfp(frf_mat(iFit,:),exp.w(iFit),Nmodes+options.nExtraModes);

            for i = 1:NHam
                for k = 1:NAccel
                    modes.fit.H(iFit,i,k) = modal_par.H(:,i+(k-1)*NHam);
                end
            end

            for j = 1:Nmodes
                [~,ii] = min(abs(modal_par.omega - setup.wEst(j)));
                modes.omega(j) = modal_par.omega(ii);
                modes.zeta(j)  = modal_par.zeta(ii);
                modal_par.omega(ii) = Inf;
                for i = 1:NHam
                    for k = 1:NAccel
                        modes.A(j,i,k) = modal_par.A(i+(k-1)*NHam,ii);
                    end
                end
            end
        otherwise
            if strcmp(options.method,'rfp')
                for i = 1:NHam
                    %node that the frequency base may be different between tests
                    for k = 1:NAccel
                        iFit = exp.w > setup.wMin & exp.w < setup.wMax;
                        modal_par = rfp(exp.H(iFit,i,k),exp.w(iFit),Nmodes+options.nExtraModes);
                        modes.fit.H(iFit,i,k) = modal_par.H;
                        
                        %now extract modal parameters
                        for j = 1:Nmodes
                            [~,ii] = min(abs(modal_par.omega - setup.wEst(j)));
                            modes.peak(i,j,k).wr = modal_par.omega(ii);
                            modes.peak(i,j,k).zr = modal_par.zeta(ii);
                            modes.peak(i,j,k).Ar = modal_par.A(ii);
                        end
                    end
                end
            else
                for i = 1:NHam
                    %node that the frequency base may be different between tests
                    for k = 1:NAccel
                        for j = 1:Nmodes
                            %find the model.peaks in the specified range
                            iBand = exp.w(2:end-1) > setup.wBand(j,1) & exp.w(2:end-1) < setup.wBand(j,2);
                            modes.peak(i,j,k) = feval(options.method,exp.w(iBand),exp.H(iBand,i,k));
                            
                            %compute fit of mode
                            modes.fit.H(iBand,i,k) = modes.peak(i,j,k).Ar ./ (modes.peak(i,j,k).wr^2 + 2*1i*modes.peak(i,j,k).zr*modes.peak(i,j,k).wr*exp.w(iBand) - exp.w(iBand).^2);
                        end
                    end
                end
            end
            
            disp('Extracting natural frequencies and modal damping...')
            modes = modal_average(modes,exp,setup);
    end
    
%     % sort modes by frequency
%     [modes.omega,iSort] = sort(modes.omega);
%     modes.zeta  = modes.zeta(iSort);
%     modes.A  = modes.A(iSort,:,:);
    
    % work out the mode shapes from the modal constants
    disp('Extracting modeshapes..')
    modes.u = modal_shapes(modes,setup.geom);
    modes.r = setup.geom.r;
    modes.n = setup.geom.n;
    modes.iBody = setup.geom.iBody;
    
    % Compute residual stiffness/mass terms if necessary
    modes.resid.K = zeros(NHam,NAccel);
    modes.resid.M = zeros(NHam,NAccel);
    for k = 1:NAccel
        for i = 1:NHam
            u = zeros(Nfreq,1);
            for j = 1:Nmodes
                u = u + modes.u(setup.geom.iHam(i),j) * modes.u(setup.geom.iAcc(k),j)./ (-exp.w.^2 + 2*1i*modes.omega(j)*modes.zeta(j)*exp.w + modes.omega(j)^2);
            end
            
            iFit = exp.w > setup.wMin & exp.w < setup.wMax;
            if options.bResidMass && options.bResidStiffness
                [modes.resid.K(i,k),modes.resid.M(i,k)] = residual_mass_and_stiffness(exp.w(iFit)',real(u(iFit)),real(exp.H(iFit,i,k)));
            elseif options.bResidStiffness
                modes.resid.K(i,k) = residual_stiffness(exp.w(iFit),real(u(iFit)),real(exp.H(iFit,i,k)));
                modes.resid.M(i,k) = Inf;
            elseif options.bResidMass
                modes.resid.K(i,k) = Inf;
                modes.resid.M(i,k) = residual_mass(exp.w(iFit),real(u(iFit)),real(exp.H(iFit,i,k)));
            else
                modes.resid.K(i,k) = Inf;
                modes.resid.M(i,k) = Inf;
            end
        end
    end
    
    save(modes_mat_file,'-struct','modes');
else
    modes = load(modes_mat_file);
end

% %plot fit at each mode
% for i = 1:NHam
%     for k = 1:NAccel
%         plot(han.exp.axMag(k),exp.w/2/pi,abs(modes.fit.H(:,i,k)),'color',options.test_col(i,:),'LineStyle','--')
%         plot(han.exp.axPh(k) ,exp.w/2/pi,angle(modes.fit.H(:,i,k))*180/pi,'color',options.test_col(i,:),'LineStyle','--')
%     end
% end


%plot the modal constants on FRF
for j = 1:Nmodes
    for i = 1:NHam
        for k = 1:NAccel
            scale = (2*1i*modes.omega(j)^2*modes.zeta(j));
            Hmax = modes.A(j,i,k)./scale;
            plot(han.exp.axMag(k),modes.omega(j)/2/pi,abs(Hmax),   'x','color',options.test_col(i,:));
            plot(han.exp.axPh(k) ,modes.omega(j)/2/pi,angle(Hmax)*180/pi,'x','color',options.test_col(i,:));
        end
    end
end

%plot natural frequencies on FRF
for k = 1:NAccel
    for j = 1:Nmodes
        han.exp.hResMag(j,k) = plot(han.exp.axMag(k),modes.omega(j)/2/pi*[1 1],ylim(han.exp.axMag(k)),'color',options.mode_col(j,:));
        han.exp.hResPh(j,k)  = plot(han.exp.axPh(k) ,modes.omega(j)/2/pi*[1 1],ylim(han.exp.axPh(k)) ,'color',options.mode_col(j,:));
    end
end


%% Compute synthesised FRFs
model_mat_file = fullfile(dataroot,'model.mat');
if ~isfile(model_mat_file) || moddate(modes_mat_file) > moddate(model_mat_file)
    model.w = exp.w;
    model.TestLabel = exp.TestLabel;
    for i = 1:NHam
        for k = 1:NAccel
            model.H(:,i,k) = - 1./(modes.resid.M(i,k)*model.w.^2) + 1/modes.resid.K(i,k);
            for j = 1:Nmodes
                model.H(:,i,k) = model.H(:,i,k) + modes.u(setup.geom.iHam(i),j) * modes.u(setup.geom.iAcc(k),j)./ (-model.w.^2 + 2*1i*modes.omega(j)*modes.zeta(j)*model.w + modes.omega(j)^2);
            end
        end
    end
    save(model_mat_file,'-struct','model')
else
    model = load(model_mat_file);
end

% %overlay model fit with experimental data
% for i = 1:NHam
%     for k = 1:NAccel
%         plot(han.exp.axMag(k),exp.w/2/pi,abs(model.H(:,i,k)),'color',options.test_col(i,:),'LineStyle','--')
%         plot(han.exp.axPh(k) ,exp.w/2/pi,angle(model.H(:,i,k))*180/pi,'color',options.test_col(i,:),'LineStyle','--')
%     end
% end

han.model = modal_plot_frfs(model,setup,options);

figure('Name',['Comparison: ' setup.Name]);

for i = 1:NHam
    for k = 1:NAccel       
        axCompare(i,k) = subplot(max(NHam,8),NAccel,(i-1)*NAccel+k);
        yyaxis left
        hold on
        if i == 1
            title(setup.AccName{k});
        end
        if k == 1
            ylabel(sprintf('%s\n%s',exp.TestLabel{i},'Mag (m/N)'));
        end
        set(axCompare(i,k),'yscale','log')
        
        if i == NHam, xlabel('f (Hz)'),  end
        plot(exp.w/2/pi,abs(exp.H(:,i,k)));
        plot(model.w/2/pi,abs(model.H(:,i,k)));
        
        yyaxis right
        phExp = angle(exp.H(:,i,k));
        phModel = unwrap(angle(model.H(:,i,k)));
        phExp = phModel + wrapToPi(phExp - phModel);
        plot(exp.w/2/pi,180/pi*phExp);
        plot(model.w/2/pi,180/pi*phModel);
        if i == NHam, xlabel('f (Hz)'),  end
        
        if k == NAccel
            ylabel('Phase (deg)')
        end
    end
end

linkaxes(axCompare(:),'x');
legend(axCompare(end),{'Exp','Model'},'AutoUpdate','off')
xlim(axCompare(1),[setup.wMin setup.wMax]/2/pi);

function [K,M] = residual_mass_and_stiffness(w,u,y)
c = [0*w+1 -1./(w.^2)]\(y-u);
K = 1/c(1);
M = 1/c(2);

function K = residual_stiffness(w,u,y)
K = 1/mean(y-u);

function M = residual_mass(w,u,y)
a = (-1./(w.^2))\(y-u);
M = -1./a;

function options = default_options(options,Ntest,Nmodes)
if ~isfield(options,'method')
    options.method = 'peak_fit';
end
if ~isfield(options,'bResidStiffness')
    options.bResidStiffness = false;
end
if ~isfield(options,'bResidMass')
    options.bResidMass = false;
end
if ~isfield(options,'nExtraModes')
    options.nExtraModes = 0;
end
if ~isfield(options,'test_col')
    options.test_col = lines(Ntest);
end

if ~isfield(options,'mode_col')
    options.mode_col = lines(Nmodes);
end