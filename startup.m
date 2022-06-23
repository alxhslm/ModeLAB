% startup script to make Matlab aware of the ModeLAB package
% Time-stamp: <2022-06-23 18:32:37 averter>

disp(['executing ModeLAB startup script...']);

me = mfilename;                                           % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));      % where am I located

addpath(mydir(1:end-1))
addpath([mydir,'Peak fitting'])
addpath([mydir,'SignalCalc'])
addpath([mydir,'Utilities'])
addpath([mydir,'Example'])

clear me mydir
