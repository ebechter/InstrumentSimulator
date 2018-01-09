function [p1] = Cheby_Fit(x,y,order,npix)
    
%     Fits Chebyshev polynomials to y as a function of x
    med = .5*npix;
    
%     bounds = [-25 25];
%     range = bounds(2)-bounds(1);
%     x0 = (x-bounds(1))/range;
    
      x_norm = x./med;
      

    
% %     normalize x
%     x_norm = (x-med)./ med;

    
    p0 = zeros( 1, order + 1 );
    p0(order) = mean( y );
    chebs = get_chebs( x_norm, order);
     
    

    
func = @(p)(errfuncCheb(p,chebs,y,order));
   
% errfuncCheb (p,chebs,y,order)
options=optimset('Display', 'off');   
    
   p1 = lsqnonlin(func,p0,[],[],options); 
   
end 
   
   
% tic
% pv_out_params = lsqnonlin(sumpv,p0,ll,ul,options);  
   
   
   
   
   
%     p1, success = scipy.optimize.leastsq(errfunc, p0, args=(chebs, y, order))
    









    