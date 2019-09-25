function frf = modal_load_frfs(dataroot)

%find the ascii files
d = dir(fullfile(dataroot,'*','Hx*.txt'));
frf_ascii_files = strcat({d(:).folder},filesep,{d(:).name})';

Data = read_signalcalc(frf_ascii_files);
NTest = size(Data,1);
NSig = size(Data,2);
switch Data(1).Type
    case 'Time'
        Fs = 1/mean(diff(Data(1).Time));
        for i = 1:NTest
            for k = 2:NSig
                [TF,Freq] = tfestimate(Data(i,1).Real,Data(i,k).Real,[],[],[],Fs);
                V(i,k-1).H = TF;
                V(i,k-1).Frequency = 2*pi*Freq;
                V(i,k-1).Real = real(TF);
                V(i,k-1).Imaginary = imag(TF);
                V(i,k-1).Magnitude = abs(TF);
                V(i,k-1).Phase = angle(TF);
                V(i,k-1).Units = Data(i,k).Units;
            end
            [V(i,:).Label] = deal(Data(i,1).Label);
        end
    case 'TF'
        V = Data;
    case 'Spectrum'
        NSig = size(Data,2);
        for i = 1:NTest
            for k = 2:NSig
                V(i,k-1).H = Data(i,k).H ./ Data(i,1).H;
                V(i,k-1).Frequency = Data(i,k).Frequency;
                V(i,k-1).Real = real(V(i,k-1).H);
                V(i,k-1).Imaginary = imag(V(i,k-1).H);
                V(i,k-1).Magnitude = abs(V(i,k-1).H);
                V(i,k-1).Phase = angle(V(i,k-1).H);
                V(i,k-1).Units = Data(i,k).Units;
            end
            [V(i,:).Label] = deal(Data(i,1).Label);
        end
end

frf.w = V(1).Frequency + eps;
Nfreq = length(frf.w);

frf.H = zeros(Nfreq,NTest,NSig);
for k = 1:NSig
    %convert acc/vel to disp
    if strncmp(V(1,k).Units,['m/s' char(178)],4)
        fprintf('Changing units of channel %d from ''m/s2'' to ''m''\n',k)
        scale = -frf.w.^2;
        units = ['m' V(1,k).Units(5:end)];
    elseif strncmp(V(1,k).Units,'m/s',3)
        fprintf('Changing units of channel %d from ''m/s'' to ''m''\n',k)
        scale = 1i*frf.w;
        units = ['m' V(1,k).Units(4:end)];
    else
        scale = 0*frf.w + 1;
        units = V(1,k).Units;
    end
    
    for i = 1:NTest
        %flip sign if hammer in the negative direction
        frf.H(:,i,k) = V(i,k).H ./ scale;
    end
end

%get legend entries
frf.TestLabel = {V(:,1).Label};