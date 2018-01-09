function ret = JointFitFunc(p,chebs,m,nx,nm)
ret = (1./m) .* Joint_Polynomial_Cheby(p,chebs,nx,nm);
end