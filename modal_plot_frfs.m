function han = modal_plot_frfs(frf,setup,options)

han.fig = figure('Name',setup.Name);
NAccel = size(frf.H,3);
Ntest = size(frf.H,2);

for k = 1:NAccel       
    han.axMag(k) = subplot(2,NAccel,0+k);
    hold on
    if k == 1
        ylabel(han.axMag(k),'Mag')
    end
    set(han.axMag(k),'yscale','log')
    title(setup.AccName{k})
    
    han.axPh(k) = subplot(2,NAccel,NAccel+k);  
    hold on
    if k == 1
        ylabel(han.axPh(k),'Ph (deg)')
    end
    xlabel('f (Hz)')
end

for i = 1:Ntest
    for k = 1:NAccel
        han.hMag(i,k) = plot(han.axMag(k),frf.w/2/pi,abs(frf.H(:,i,k)),'-','color',options.test_col(i,:));
        han.hPh(i,k)  = plot(han.axPh(k) ,frf.w/2/pi,angle(frf.H(:,i,k))*180/pi ,'color',options.test_col(i,:));
        
        if setup.bDrivePt(i,k)
            set(han.hMag(i,k),'LineWidth',2)
            set(han.hPh(i,k),'LineWidth',2)
        end
    end
end

linkaxes([han.axMag(:);han.axPh(:)],'x');
legend(han.axMag(end),han.hMag(:,end),frf.TestLabel,'AutoUpdate','off')
xlim(han.axMag(1),[setup.wMin setup.wMax]/2/pi);
