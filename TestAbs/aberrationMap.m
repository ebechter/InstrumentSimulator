% interpolate per order per aberration. 

nOrders = 10;
nAberrations = 5;

% abberArray = (9x10x36x3)



for ii = 1:nOrders
    
    aberrSamples = abberArray(:,ii);
    wavelengthSamples = wavelengthSampleArray(:,ii);
    wavelengthVector = wavelengthArray(:,ii);
    
    for jj = 1:nAberrations
        a{ii}(:,jj) = interp1(wavelengthSamples,aberrSamples,wavelengthVector);
    end
end