function ret_val = Cheby_eval(p,x,npix)
% 
% evaluates Chebyshev polynomial fit at x given best-fit parameters p
% 
med = .5*npix;
x_norm = x./med;

% x_norm = (x-med)./ med;
order = length(p) - 1;
ret_val = 0.0;
for i = 0:order %in range(order + 1):
    ret_val = ret_val+ p(order+1 - i)* chebyshevT(i,x_norm);   % chebyscipy.special.chebyt(i)(x_norm)
end
end