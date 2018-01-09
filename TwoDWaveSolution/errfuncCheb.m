function out = errfuncCheb (p,chebs,y,order) 
    out = fitfuncCheb(p,chebs,order)-y; 
end