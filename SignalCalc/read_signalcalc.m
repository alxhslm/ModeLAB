function V = read_signalcalc(file)
% READ_SIGNALCALC_FRF reads in SignalCalc ASCII data files
%
% Inputs:
%       file: name (or full path) of file to be read
%
% Outputs:
%       V: structure array containing the data, where the length is equal to
%          the number of channels. It will contain the following fields:
%           - Frequency or Time
%           - Real
%           - [Imaginary]
%           - [Magnitude]
%           - [Phase]

if ~iscell(file)
    [~,label] = fileparts(file);
    fid = fopen(file);

    %read first line to determine what format we have
    hdr = fgetl(fid);
    if strncmp(hdr,'Data written for',16)
        %multiple channel format
        [sig,chan,units,data] = multi_chan(fid);
    else
        %single channel format
        [sig,chan,units,data] = single_chan(fid,label(1));
    end

    switch label(1)
        case {'X','W'}
            type = 'Time';
        case 'H'
            type = 'TF';
        case 'S'
            type = 'Spectrum';
        otherwise
            type = 'Unknown';
    end
    
    %sometimes some channels share units
    NSig = length(sig);
    unitsX = units(1);
    unitsY = reshape(units(2:end),[],NSig);
    NChan = size(unitsY,1);
    for k = 1:NSig
        for l = 1:NChan
            if strcmp(unitsY(l,k),'')
                unitsY(l,k) = unitsY(1,k);
            end
        end
    end
    units = [unitsX unitsY(:)'];

    %convert to SI units
    for i = 1:length(chan)
        switch units{i}
            case 'sec'
                scale = 1;
                units{i} = 's';
            case 'Hz'
                scale = 2*pi;
                units{i} = 'rads';
            case 'deg'
                scale = pi/180;
                units{i} = 'rad';
            otherwise
                scale = 1;
        end
      data(:,i) = scale*data(:,i);
    end

    %extract common X vector
    chanX = chan{1};
    dataX = data(:,1);
    unitsX = units{1};

    %extract Y data
    chanY = reshape(chan(2:end),NChan,NSig);
    dataY = data(:,2:end);
    unitsY = reshape(units(2:end),NChan,NSig);
    
    if strcmp(type,'Time') && NSig > 2
        %signalcalc exports time series with multiple channels wrongly
        %fix this issue here
        iTrim = find(~isnan(dataY(:,3)),1)-1;
        dataY(:,3:end) = [dataY((iTrim+1):end,3:end); dataY(end-iTrim+1:end,1:(NSig-2))];
        dataY = dataY(1:end-iTrim,:);
        dataX = dataX(1:end-iTrim,:);
    end
    
    for k = 1:NSig
        V(k).(chanX) = dataX;
        for l = 1:NChan
            V(k).(chanY{l,k}) = dataY(:,(k-1)*NChan+l);
            V(k).Units = unitsY{l,k};
        end
        V(k).Name = sig{k};
        V(k).Label = label;
        V(k).Type = type;
    end
    
    %finally convert mag/ph -> re/im and vice versa
    if ~isfield(V(1),'Magnitude') && (isfield(V(1),'Real') && isfield(V(1),'Imaginary'))
        for k = 1:NSig
            [V(k).Phase,V(k).Magnitude] = cart2pol(V(k).Real,V(k).Imaginary);
        end
    elseif ~isfield(V(1),'Real') && (isfield(V(1),'Magnitude') && isfield(V(1),'Phase'))
        for k = 1:NSig
            [V(k).Real,V(k).Imaginary] = pol2cart(V(k).Phase,V(k).Magnitude);
        end
    end

    %create complex channel
    if ~isfield(V,'Imaginary')
        for k = 1:NSig
            V(k).Imaginary = 0*V(k).Real;
        end
    end
    
    for k = 1:NSig
        V(k).H = V(k).Real + 1i*V(k).Imaginary;
    end

    if isfield(V(1),'Phase')
        for k = 1:NSig
            V(k).Phase = unwrap(V(k).Phase);
        end
    end
    
else
    for i = 1:length(file)
        V(i,:) = read_signalcalc(file{i});
    end
end

function [sig,chan,units,data] = single_chan(fid,label)
%extract complete header
i = 1;
hdr{1} = fgetl(fid);
while ~strcmp(hdr{i},'Units:')
    i = i + 1;
    hdr{i} = fgetl(fid);
end
%one more for units
hdr{end+1} = fgetl(fid);

%no signal name for single channel files
sig = {[label '1']};
NSig = 1;

%extract channel names
hdrChan = hdr{end-3};
chan = strtrim(strsplit(hdrChan,'\t')); 
NChan = length(chan)-1;

%extract units
hdrUnits = hdr{end};
units = strtrim(strsplit(hdrUnits,'\t'));
if length(units)<NChan+1
    units{NChan+1} = '';
end

%and finally the actual data
data  = textscan(fid,repmat('%f' ,1,NChan*NSig+1),'HeaderLines',1,'CollectOutput',1); data = data{1};
fclose(fid);

%rather stupidly, for files with a single channel, SignalCalc writes ratio
%units as (m/s)(N) rather than m/s/N like in files with multiple channels
unitsX = units(1);
unitsY = units(2:end); 
for i = 1:length(unitsY)
    tmp = regexp(unitsY{i},'(?<=\()(.*?)(?=\))','match');
    if length(tmp)>1
        unitsY{i} = [tmp{1} '/' tmp{2}];
    elseif isempty(tmp)
        unitsY{i} = '';
    else
        unitsY{i} = tmp{1};
    end
end
units = [unitsX unitsY];

function [sig,chan,units,data] = multi_chan(fid)
%extract complete header
i = 1;
hdr{1} = fgetl(fid);
while ~strcmp(hdr{i},'Units:')
    i = i + 1;
    hdr{i} = fgetl(fid);
end
%one more for units
hdr{end+1} = fgetl(fid);

%list of signals - eg H1,2  H1,3 etc
hdrSig = hdr{end-5};
sig = strtrim(strsplit(hdrSig,'\t')); sig = sig(2:end);  sig = sig(~strcmp(sig,'')); 
NSig = length(sig);

%extract channel names - eg Real Imaginary Magnitude
%likely >1 per signal
hdrChan = hdr{end-4};
chan = strtrim(strsplit(hdrChan,'\t'));
NChan = (length(chan)-1)/NSig;

%extract units eg m/s2 deg etc
hdrUnits = hdr{end};
units = strtrim(strsplit(hdrUnits,'\t'));
if length(units)<NChan*NSig+1
    units{NChan*NSig+1} = '';
end

%and finally the actual data
data  = textscan(fid,repmat('%f' ,1,NChan*NSig+1),'HeaderLines',1,'CollectOutput',1); data = data{1};
fclose(fid);

% rather stupidly, for files with multiple channels, SignalCalc calls it 
% Imag instead of Imaginary 
[chan{strcmp(chan,'Imag')}]= deal('Imaginary');
[chan{strcmp(chan,'Mag')}]= deal('Magnitude');