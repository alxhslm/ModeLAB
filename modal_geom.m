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