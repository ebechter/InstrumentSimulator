function [m0_ret] = find_m(lamc)


% lamc =   % column vector
N=200;

m0 = 1:N;

j = flipud((1:36)');


% assert(size(j)==size(lamc));

m0lj = lamc*m0;

jla = repmat(j.*lamc,1,length(m0));

[~,ind]=min(mean(abs(diff(m0lj+jla)),1));

m0_ret = ind;
