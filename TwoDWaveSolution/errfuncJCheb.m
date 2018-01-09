function out = errfuncJCheb (p,chebs,y,m,nx,nm) 
    out = JointFitFunc(p,chebs,m,nx,nm)-y; 
end