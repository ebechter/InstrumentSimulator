function out = fitfuncCheb(p,chebs,order)
ret_val = 0.0;
for ii = 0:order
    ret_val = ret_val+ p(order+1-ii)*chebs{ii+1};
    out =ret_val;
end