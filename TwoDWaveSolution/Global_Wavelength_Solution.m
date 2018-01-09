function [p1, ret, chebs] = Global_Wavelength_Solution(pix_centers, wavelengths, orders, p0,...
    order0, ntotal, npix, Inverse,nx,nm)


chebs = Calculate_chebs(pix_centers, orders+order0, Inverse,order0,ntotal,npix,nx,nm);


%     p0 = zeros( 1, order + 1 );
%     p0(order) = mean( y );

    
func = @(p)(errfuncJCheb(p,chebs,wavelengths,orders+order0,nx,nm));


% errfuncCheb (p,chebs,y,order)
options=optimset('Display', 'off');   
p1 = lsqnonlin(func,p0,[],[],options); 


% p1, success =  scipy.optimize.leastsq(errfunc_cheb, p0, args=(chebs, wavelengths, orders + order0, Wgt))
% residuals    = errfunc_cheb_nw(p1, chebs, wavelengths, orders + order0)
        
%         residuals_ms = 299792458.0 * residuals / wavelengths;
%         rms_ms     = np.sqrt( np.var( residuals_ms ) );
%         
        
%         figure
%         hold on
%     
%         for ii = 1:39
%         plot
%         erf = ret-wavelengths;
      

        
ret = JointFitFunc(p1,chebs,orders+order0,nx,nm);


        
end
    

