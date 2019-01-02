
FRF.mat contain:
    FRF(:,1) = Frequency range vector (rad/sec).
    FRF(:,2) = The simulated Frequency Response Function using as input=1 and output=3
             (see: MATLAB Central>File Exchange>Controls and Systems Modeling>Mechanical Modeling>Linear forced system with viscous damping).
   

FRF_noise.mat contain:
    FRF_noise(:,1) = Frequency range vector (rad/sec).
    FRF_noise(:,2) = A noisy Frequency Response Function.
    noise was generated as:
        >> IRF=ifft(FRF,512);
        >> max_irf=max(abs(IRF));
        >> noise=0.005*max_irf*randn(size(IRF));
        >> sprintf('SNR = %0.5g [dB].',20*log10(std(IRF)/std(noise))),
        >> IRF=IRF+noise;
        >> FRF_noise=fft(IRF,512); %plot(w,20*real(log10(FRF_noise)),'b'),

 To use the rfp.m function, you must do the following steps:
 1) Choose the frequency range of the FRF. This vector is called 'omega' and the respective values of the FRF is called 'rec'.
 2) Choose the number of modes in the range (in this case N=3).
 3) Run rfp.m using "rec, omega and N".
 You will obtain the estimated FRF (alpha vector). 
 Plot the estimated FRF (alpha) with the simulated FRF (FRF).
      
 i.e.:  >> load FRF_noise.mat
        >> plot(FRF_noise(:,1),20*log10(abs(FRF_noise(:,2))),'b'), hold on,
        >> xlabel('Frequency (rad/sec)'),ylabel('Magnitude (dB)'),
        >> omega=FRF_noise(50:450,1); 
        >> rec=FRF_noise(50:450,2);
        >> N=3; 
        >> [alpha,par]=rfp(rec,omega,N);
        >> fn=par(:,1) %natural frequencies
        >> xi=par(:,2) %damping ratios
        >> C=[par(:,3),par(:,4)] %modal constant (magnitud,phase)
        >> plot(omega,20*log10(abs(alpha)),'r'), hold off,
