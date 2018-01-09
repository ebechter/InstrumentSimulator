clear
wfe = [0,0,0,0,0,0,0,0]; % Pupil Zernike amplitudes piston, xtilt, ytilt, etc...
rms_scale = [1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8),1/sqrt(8)];
adc = 0; %microns is max ADC dispersion (This is translated as an offset fiber pos rather than beam position for computational ease)
fiberpos = [0,0,0]; % global position offset in microns (x,y,z)
dof = 0;
amp = 0.5;

figure
hold on
for jj = 2
    wfe = [];
    nc = [];
    c = [];
for ii = 1:length(amp)
    wfe = [0,0,0,0,0,0,0,0];
    fiberpos = [0,0,0];
    R = FiberCoupling(wfe,adc,fiberpos,dof);
    c(:,ii) = R.Rho
    rms_t = wfe.*rms_scale;
end
plot(c)
end

% 
% PSF(:,:,1) = R.FPgridx(:,:,1);
% PSF(:,:,2) = R.FPgridy(:,:,1);
% PSF(:,:,3) = abs(R.PSF(:,:,1)).^2;

PSF(:,:,1) = R.FPgridx(:,:,1);
PSF(:,:,2) = R.FPgridy(:,:,1);
PSF(:,:,3) = abs(R.FiberPSF(:,:,1)).^2;