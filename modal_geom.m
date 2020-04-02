function geom = modal_geom(setup)
%identify how many unique points we have
rLoc = [setup.rHam;
        setup.rAcc];
    
nLoc = [setup.nHam;
        setup.nAcc]; 

A = [rLoc nLoc];
[C,iA,iC] = unique(A,'rows','stable');

%store unique set of points
geom.r = C(:,1:3);
geom.n = C(:,4:6);

NHam = size(setup.rHam,1);

%find which unique point each hammer/accel location corresponds to
geom.iHam = iC(1:NHam);
geom.iAcc = iC(NHam+1:end);

%and which unique points come from hammer/accel
geom.iFromHam = iA(iA <= NHam);
geom.iFromAcc = iA(iA > NHam) - NHam;

%finally, on which body is each unique point
geom.iBody = [setup.iBodyHam(geom.iFromHam); setup.iBodyAcc(geom.iFromAcc)];

%find driving point responses (where iHam == iAcc)
geom.bDrivePt = (geom.iHam - geom.iAcc') == 0;

for i = 1:size(setup.nAcc,1)
   geom.bHamAccParallel(:,i) =  setup.nHam * setup.nAcc(i,:)' == 1;
   geom.bHamAccSameBody(:,i) = setup.iBodyHam == setup.iBodyAcc(i);
end

%check if each test is in the direction of each mode
if isfield(setup.modes,'nMode')
    for k = 1:setup.modes.Nmodes
        geom.bModeAcc(:,k) = setup.nAcc * setup.modes.nMode(k,:)' ~= 0;
        geom.bModeHam(:,k) = setup.nHam * setup.modes.nMode(k,:)' ~= 0;
    end
else
    geom.bModeAcc = true(size(setup.rAcc,1),setup.modes.Nmodes);
    geom.bModeHam = true(size(setup.rHam,1),setup.modes.Nmodes);
end
