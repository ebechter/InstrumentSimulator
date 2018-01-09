function [coefs] =Calculate_chebs(x,m, order0, ntotal,npix,Inverse,nx,nm)
% #normalize x

med = 0.5*npix;
x_norm = x./ med;



if Inverse == 1
%     # invert and normalize m
    im = 1/m;
    im_max = 1./order0;
    im_min = 1./order0+ntotal;
else
    im = m;
    im_max = order0+ntotal;
    im_min = order0;
end
    delta = 0.5*(im_max - im_min);
    m_norm = (im-im_min-delta)./delta;


    coefs = [];

    for ii = 1:nx
        
        coefs{ii} = chebyshevT(ii,x_norm);
    end
    
    for ii = 1:nm
        
        coefs{nx+ii} = chebyshevT(ii,m_norm);
    end 
end
%             """
%             u = sp.chebyt(1)(x_norm)
%             u2 = sp.chebyt(2)(x_norm)
%             u3 = sp.chebyt(3)(x_norm)
%             u4 = sp.chebyt(4)(x_norm)
%             v  = sp.chebyt(1)(m_norm)
%             v2 = sp.chebyt(2)(m_norm)
%             v3 = sp.chebyt(3)(m_norm)
%             v4 = sp.chebyt(4)(m_norm)
%             v5 = sp.chebyt(5)(m_norm)
%             v6 = sp.chebyt(6)(m_norm)
%             return u,u2,u3,u4,v,v2,v3,v4,v5,v6
%             """
