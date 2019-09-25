function setup = modal_setup(setup_csv_file)
P = read_setup_csv(setup_csv_file);

%% Geometry
setup.nAcc = [P.axAccel;
             P.ayAccel;
             P.azAccel]';
setup.nHam = [P.FxHammer;
             P.FyHammer;
             P.FzHammer]';
setup.rHam = [P.xHammer;
             P.yHammer;
             P.zHammer]';
setup.rAcc = [P.xAccel;
             P.yAccel;
             P.zAccel]';
         
NHam = size(setup.rHam,1);
NAcc = size(setup.rAcc,1);

%flip sign if accelerometer in the negative direction
sHam = (1-2*any(setup.nHam < 0,2));
sAcc = (1-2*any(setup.nAcc < 0,2));
setup.sTest = sHam*sAcc';

setup.nHam = setup.nHam .* sHam;
setup.nAcc = setup.nAcc .* sAcc;

setup.iBodyHam = P.iHammer';
setup.iBodyAcc = P.iAccel';

%identify how many unique points we have
rLoc = [setup.rHam;
        setup.rAcc];
    
nLoc = [setup.nHam;
        setup.nAcc]; 

A = [rLoc nLoc];
[~,~,iC] = unique(A,'rows','stable');

%find index of each hammer/accel location
iHam = iC(1:NHam);
iAcc = iC(NHam+1:end);

%find driving point responses (where iHam == iAcc)
setup.bDrivePt = false(NHam,NAcc);
for k = 1:NAcc
    ii = find(iHam == iAcc(k));
    if ~isempty(ii)
        setup.bDrivePt(ii,k) = true;
    end
end

%% Frequency information
setup.wBand = [P.fL' P.fH']*2*pi;
setup.wEst = mean(setup.wBand,2);

if isfield(P,'fMin')
    setup.wMin = P.fMin*2*pi;
else
    setup.wMin = 0;
end

if isfield(P,'fMax')
    setup.wMax = P.fMax*2*pi;
else
    setup.wMax = Inf;
end

%% Accelerometer labels
if ~isfield(P,'AccelName')
    for k = 1:length(P.axAccel)
        if P.axAccel(k) ~= 0 && P.azAccel(k) == 0
            setup.AccName{k} = 'Horizontal';
        elseif P.ayAccel(k) ~= 0 && P.azAccel(k) == 0
            setup.AccName{k} = 'Vertical';
        elseif P.azAccel(k) ~= 0 && P.axAccel(k) == 0 && P.ayAccel(k) == 0
            setup.AccName{k} = 'Axial';
        else
            setup.AccName{k} = 'Unknown';
        end
    end
else
    setup.AccName = P.AccelName;
end
