options.method = 'peak_fit';
options.bPlotModel = false;
options.bPlotComparison = false;
    
modes = modal_analysis('./Example/SignalCalc',options);

figure
for i = 1:4
    subplot(4,1,i)
    plot(modes.r(:,3), modes.u(:,i))
    title(sprintf('%0.2f Hz', modes.omega(i)/2/pi))
end