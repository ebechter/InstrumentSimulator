exponent = ((x).^2 + (y).^2)./(gaussfit.c1^2);
amplitude = gaussfit.a1;  
% The above is very much different than Alan's "1./2*pi*sigma^2"
% which is the same as pi*Sigma^2 / 2.
val       = amplitude  * exp(-exponent);




