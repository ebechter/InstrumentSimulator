function ret_val = Joint_Polynomial_Cheby(p,chebs,nx,nm)
%
%     """
%     Evaluates a *product* of Chebyshev polynomials in x and m
%     Polynomial is tailored.
%     """

xvec = chebs(1:nx);  % might need to change this indexing... 
mvec = chebs(nx+1:end);
ret_val = p(1);
k=2;
for ii = 1:nx
    ret_val = ret_val + p(k).*xvec{ii};
    k= k+1;
end
for ii = 1:nm
    ret_val = ret_val + p(k).*mvec{ii};
    k= k+1;
end

if nx >= nm
    for ii = 1:nx
        for jj = 1 : min(nx-ii,nm)
            
            ret_val = ret_val+ p(k).*xvec{ii}.*mvec{jj};
            k=k+1;
        end
    end
            
else
    for jj =1:nm
        for ii = 1:min(nm-jj-1,nx)
            ret_val = ret_val + p(k)*xvec{ii}.*mvec{jj};
            k= k+1;
        end
    end
    
    
    
end