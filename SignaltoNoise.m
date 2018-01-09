function[out]=SignaltoNoise(bandpass,flux_raw)
ii = 1;
integration_time = 30*60; %measured in seconds
Instrument = 1;% 1 = iLocater 2 = Hamilton

% fine = [0.0000000001:0.000001:0.005];
% coarse = [0.005+0.000001:0.001:0.01];
% range1 = [fine coarse];

range1 = 0.0000000001:0.000001:0.005;
range2 = [1E-7,0.05]

tput_analysis = 0;
% range2 =0.015;
etalon = 5000;

for tput =  range2
    flux = flux_raw*tput;
    %----------Optical Count Conversions----------%
    Counts_bin = flux; % multiply Count Desnity by step to get count per wavelength bin
    %Counts bin is measured in photons/sec (for a single sample)
    cnts_full = Counts_bin*integration_time; % accounts for integration time
    
    %----------Adjust for iLocater resolution and bandpass----------%    
    %%Check by cutting and comparing area%%
    % wavelength_limits = bandpass<=ru & bandpass>=rl; % Instrument band
    % cnts = cnts_full(wavelength_limits);% cut the flux value at the band
    % bandpass_asd = bandpass(wavelength_limits);
    % area1 = sum(cnts);
    %%end of check%%
    if Instrument == 1
        FWHM = 3;
        R = 180E3;
%         band_cut = 0.971+0.050:(1.14/(FWHM*R)):1.270-0.050;
        band_cut = 0.971:(1.14/(FWHM*R)):1.270-0.020; %20nm band lost        
    elseif Instrument == 3
        FWHM = 3;
        R = 180E3;
        band_cut = 0.971:(1.14/(FWHM*R)):1.270;    
    elseif Instrument == 2
        FWHM = 1.9355;
        R = 62E3;
        band_cut = 0.480:(0.550/(FWHM*R)):0.620;
    else
        disp('Pick a valid instrument')
    end
    
    [bandpass,B,C]=unique(bandpass);
    cnts_full = cnts_full(B);
    
    cnts_cut =interp1(bandpass,cnts_full,band_cut);
    
    %%area comapare%%
    % area2 = sum(cnts_cut);
    % area1/area2
    %%
    
    %----------PSF creation and convolution----------%
    x = (0:10);
    mu = 5;
    %FWHM = 2ln2*(Sigma)
    Sigma = FWHM/(2*sqrt(2*log(2)));
    [psf] = normpdf(x,mu,Sigma);
    N = conv(cnts_cut,psf,'Same');
    % N = cnts_cut; %if no convolustion
    c = 299792458;
    
    %----------Calculate Signal to noise----------%
    
    if Instrument ==1
        telescopes = 1/sqrt(2);
    else
        telescopes=1;
    end
    
    ave_count = (N(2:end)+N(1:end-1))/2;
    ave_wave = (band_cut(2:end)+band_cut(1:end-1))/2;
    e = sqrt(ave_count)./ave_count; %1/sqrt(2) -->2 telescopes
    
    %----------Continuum Normalize----------%
    spectrum(:,1) = band_cut;
    spectrum(:,2) = N;
    [spectrum] = PiecewiseDetrend(spectrum);
    band_cut=spectrum(:,1)';
    N=spectrum(:,2)';
    
    if Instrument == 1
    rl = 1.110;
    ru = 1.160;
    wavelength_limits = band_cut<=ru & band_cut>=rl; % Instrument band
    N(wavelength_limits)=1;
    end
    
    dI=diff(N);
%     cnts = cnts_full(wavelength_limits);% cut the flux value at the band
%     bandpass_asd = bandpass(wavelength_limits);
    
    %Calculate derivative
    dl = diff(band_cut);
    dv = dl*c./ave_wave;
    di_dv=dI./dv;
    
    tel_fraction = 1-mean([0.55,0.20]);
%     di_dv = di_dv*tel_fraction; %telluric modification
    
    A = (di_dv./(e)).^2;
    sig = (1/sqrt(sum(A)))*telescopes;
    SN = mean(1./e);

    temp(ii,1)=tput;
    temp(ii,2)=sig;
    temp(ii,3)=SN;

    ii = ii+1;
end

if tput_analysis == 1
%fine gridding
f_tput(:,1) = 0.0001:1E-6:1;
f_sigma(:,1) = interp1(temp(:,1),temp(:,2),f_tput,'linear','extrap');
% ind = find(0.29<f_sigma<0.31 ,1,'first');
% ind = find(0.59999< f_sigma <0.60001 ,1,'first');
target= 0.6;
[num ind] = (min(abs(f_sigma-target)));

%f_tput(:,1) = 0.001:0.001:1;
f_SN(:,1) = interp1(temp(:,1),temp(:,3),f_tput,'linear','extrap');
% f_SN(:,1) = interp1(temp(:,2),temp(:,3),f_sigma,'linear','extrap');
% ind2 = find(0.49< f_sigma <0.51 ,1,'first');

out(1,1)=f_sigma(ind,1);
out(1,2)=f_tput(ind,1);
out(1,3)=f_SN(ind,1);
else
out(:,1)=temp(:,1);
out(:,2)=temp(:,2);
out(:,3)=temp(:,3);    
end

end



