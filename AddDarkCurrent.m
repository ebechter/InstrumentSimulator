function [frame, sn] = AddDarkCurrent(dc_s,intTime,signalarray)

% Add HxRG dark current

nx = size(signalarray,1);
ny = size(signalarray,2);

dark_frame = poissrnd(dc_s*intTime,nx,ny);

frame = signalarray + dark_frame;

sn = max(max(signalarray))/max(max(dark_frame));

end