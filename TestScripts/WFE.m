clear
wfe = [0,0,0,0,0,0,0,0,0,0,0]; % Pupil Zernike amplitudes piston, xtilt, ytilt, etc...
rms_scale = [1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8),1/sqrt(8)];
adc = 0; %microns is max ADC dispersion (This is translated as an offset fiber pos rather than beam position for computational ease)
fiberpos = [0,0,0]; % global position offset in microns (x,y,z)
dof = 0;
amp = 0:0.1:1;

figure
hold on
for jj = 2:8
    wfe = [];
    nc = [];
    c = [];
for ii = 1:length(amp)
    wfe(1,jj) = amp(ii);
    R = FiberCoupling(wfe,adc,fiberpos,dof);
    c(:,ii) = R.Rho;
    rms(ii) = amp(ii)*rms_scale(1,jj);
end
nc = c./c(1,1);
jj
plot(rms,nc)
end