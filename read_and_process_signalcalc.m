function [V,leg] = read_and_process_signalcalc(test)

if ~iscell(test)
    test = {test};
end

Ntest = length(test);

disp('Openings FRFs...')
for i = 1:Ntest
    disp(test{i})
    V(i,:) = read_signalcalc(test{i});
    [~,leg{i},~] = fileparts(test{i});
end
% figure
% ax(1) = subplot(2,1,1);
% hold on
% ax(2) = subplot(2,1,2);
% hold on
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