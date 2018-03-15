classdef Simulation
    properties
        scale
        
    end
    
    methods
        function obj = Simulation(scale)
        
            if nargin == 0 
                obj.scale = 1;
            else
                obj.scale = scale;
            end
        end
        
    end
    
    methods(Static)
        
        function [sampledstar] = addStar(R,pixSamp,scale,wavelength,counts,rv)
            
            
            % Define wavelength step (input in nm, output in microns) using maximum resolving power (R),
            % pixSamp, and scale factor
            dlam = 1000/R/pixSamp;
            step = dlam/scale;
            interplambda = wavelength(1):step:wavelength(end);
            interpflux = interp1(wavelength, counts, interplambda,'linear');
            
            dsWavelength = Star.dopplerShift(interplambda,rv); %Shift Wavelength and assign DsWavelength
            counts = Star.energy2Counts(dsWavelength,interpflux);
            sampledstar = [dsWavelength;counts]';
            
            
        end
        
        function [counts] = addCollectingArea(countsPerArea,dTel,blockFrac)
            
            collecting_area = pi*(dTel/2)^2 - pi*(blockFrac*dTel/2)^2;
            counts = countsPerArea * collecting_area;
            
        end
        
        function [y_grid] = resampleToGrid(x,y,x_grid)
            % this function cuts a new x vector and the max of the original and
            % extrapolates the short end of the new x vector to the short end of the
            % original x vector. The new y vector for use is generated by interpolating the
            % new_x,new_y onto the original x vector.
            
            extrap_scalar = min(y);
            ub = min(x(end,1),x_grid(end,1));
            ind = find(x<=ub,1,'last');
            x_safe = x(1:ind);
            y_safe = y(1:ind);
            
            y_grid = interp1(x_safe(:,1),y_safe(:,1),x_grid,'linear',extrap_scalar);
            
        end
        
        function [addtell] = addAtmosphere(wavelength,counts, telluric, skyback)
            
            skygrid = Simulation.resampleToGrid(skyback(:,1),skyback(:,2),wavelength);
            
            addsky = counts + skygrid;
            
            tellgrid = Simulation.resampleToGrid(telluric(:,1),telluric(:,2),wavelength);
            
            addtell = addsky.*tellgrid;
            
        end
        
        function [starThroughput] = combineImagerThroughput(star_components)
            comTput{1} = star_components(1).finalThroughput;
            
            for ii = 2:size(star_components,2) %loop over cellarry size (i.e. surfaces)
                [comTput{ii}(:,1),comTput{ii}(:,2)] = Instrument.multiply_curves...
                    (comTput{ii-1}(:,1),comTput{ii-1}(:,2),...
                    star_components(ii).finalThroughput(:,1),star_components(ii).finalThroughput(:,2));
            end
            
            starThroughput = comTput{end};
        end
        
        function [output] = Xdisperse(inputSpectrum,nOrders,waveCoeff)
            
            
            temp{1} = repmat(inputSpectrum(:,1)',nOrders,1)';
            temp{2} = repmat(inputSpectrum(:,2)',nOrders,1)';
            
            
            [outputWavelength, outputOrders] = Simulation.trimOrders(temp{1}, temp{2}, waveCoeff,nOrders);
            
            
            output{1} = outputWavelength;
            output{2} = outputOrders;
            
        end
        
        function [OutWavelength, OutSpectrum] = trimOrders (wavelength,spectrum,wave_coeff,nOrders)
            
            ord = size(spectrum, 2); % what is the largert order (39)
            buffer_edge = 25; %set detector edge (buffered by 2mm)
            
            ind_lim1 = find(wavelength(:,ord) <= polyval(wave_coeff(ord,:,1),-buffer_edge),1,'last');
            ind_lim2 = find(wavelength(:,ord) >= polyval(wave_coeff(ord,:,1), buffer_edge),1,'first');
            max_length = (ind_lim2-ind_lim1); % find the length of the largest order with buffer
            
            for ii = 1:nOrders
                
                ind1 = find(wavelength(:,ii) <= polyval(wave_coeff(ii,:,1), -buffer_edge),1,'last');
                ind2 = max_length + ind1; % ind2 is forced to match largest order size
                %                 disp(ind2-ind1) %check all orders are the same size
                OutWavelength(:,ii) = wavelength(ind1:ind2,ii);
                OutSpectrum(:,ii) = spectrum(ind1:ind2,ii);
                %                 disp(ii)
            end
            
        end
        
        function [output] = addSpecThroughput(spectral_cell,finalThroughput,nOrders)
            
            for ii = 1:nOrders
                new_y = Simulation.resampleToGrid(finalThroughput{2}(:,ii)*1e-3,finalThroughput{1}(:,ii),spectral_cell{1}(:,ii));
                output{1}(:,ii) = spectral_cell{1}(:,ii);
                output{2}(:,ii) = new_y .* spectral_cell{2}(:,ii);
                
            end
            
        end
        
        function [PSF,center] = MakePSF(scale,pix_samp,vert_samp)
            
            % dispersion direction
            fwhm = pix_samp*scale;
            sigmax=fwhm/(2*sqrt(2*log(2)));
            
            % cross dispersion direction
            fwhmy = vert_samp*scale;
            sigmay=fwhmy/(2*sqrt(2*log(2)));
            
            % trim according to bigger dimension
            Nsigmas = 9;
            val = Nsigmas * sigmay;
            %
            
            %             clip = 3/(2*sqrt(2*log(2)));
            %             val = 25*clip;
            
            oddGrid=2*round(val/2)-1; % not super important but lets go with an odd grid size every time
            
            [MatX,MatY]=meshgrid(1:1:oddGrid,1:1:oddGrid); % make the mesh grid
            
            center = [round(size(MatX,2)/2),round(size(MatY,1)/2)]; % center x and y in the (likely) non-integer center of grid ;
            % this is the most important line - the center of the PSF must fall in a pixel center
            % subtract the center in the x direction from the full convolution to get a "central" conv. result
            %....maybe not certain.
            
            PSF=circ_gauss(MatX,MatY,[sigmax,sigmay],center);
            
            
        end
        
        function [trim, wavelength] = ConvolveOrder(wavelength,spectrum,wave_coeff,scale)
            
            % trim each order beyond the edge of the detector
            
            ind1 = find(wavelength <= polyval(wave_coeff,-22.48),1,'last');
            ind2 = find(wavelength >= polyval(wave_coeff, 22.48),1,'first');
            
            wavelength = wavelength(ind1:ind2)';
            spectrum = spectrum(ind1:ind2)';
            
            % sampled at the high end. 3 pixels at red, smooth function to
            % 6 pixels at blue.
            
            pix_samp = linspace(6,3,length(wavelength));
            
            vert_samp = median(pix_samp);
            pix_samp = 0.8*pix_samp;
            % Custom convolution
            
            % Do the first loop iteration outside the loop. Need to
            % calculate dim first
            ii = 1;
            
            [PSF,center] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);
            dim = size(PSF,1);
            rectangle = zeros(dim,length(wavelength)+dim-1);
            rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
            
            for ii = 2:length(wavelength)
                
                [PSF,~] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);
                
                %                 Conv1 = conv2(spectrum(ii),PSF,'full');
                
                rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
                %                 rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+Conv1;
                
                
            end
            
            
            trim = rectangle(:,center(1):end-(size(PSF,2)-center(1)));
            
        end
        
        function detector = CliptoDetector(order,wavelength,order_coeff,wave_coeff,cheby,p1,specorder)
            
            
            % convert each xypair into a set of 4 verticies for the pixel. Depends on PSF scaling etc..
            % from PSF convolution - know what one square "sample" is in mm. Does that still work after convolution?
            %% XY location pair for each point on wavelength vector
            xy_pair = zeros([size(order),2]);
            
            [~,beta] = max(order);
            center = mode(beta);
            
            % Convert lambda to X based on Zemax order fit
            % new x  =  interpolate (lam,x,lam')
            
            %             xy_pair(center,:,1) = interp1( polyval(wave_coeff,((1:0.1:4096)-2048)/100)  , ((1:0.1:4096)-2048)/100   ,wavelength,'linear','extrap');
            
            %             xy_pair(center,:,1) = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,wavelength,'linear','extrap');
            
            % Chebyshev 2D polynomial wavelength interpolation.
            specorder = 40-specorder;
            if cheby == 1
                Inverse =1;
                nx = 7;
                nm = 7;
                ntotal = 39;
                npix = 4096;
                order0 = 113;
                finex=-25.48:0.01:25.48;
                fineorders = specorder*ones(size(finex)); % what spectral order are we at?
                finechebs = Calculate_chebs(finex, fineorders+order0,Inverse,order0,ntotal,npix,nx,nm); 
                finewavelengths = JointFitFunc(p1,finechebs,fineorders+order0,nx,nm); 
                
                xy_pair(center,:,1) = interp1( finewavelengths,finex,wavelength,'linear','extrap');
                
            else
                % using 4th order, 1D polynomial needs interpolation. Inversion
                % is not possible.
                xy_pair(center,:,1) = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,wavelength,'linear','extrap');
            end
            
            % Inverting the 2nd order 1D polynomial
            %             a = wave_coeff(1);
            %             b = wave_coeff(2);
            %             c = wave_coeff(3);
            %
            %             xy_pair(center,:,1) = (-b + sqrt(b^2-4*a*(c-wavelength)))/(2*a);
            
            
            % fprintf('Counts before detector: %.5f', sum(sum(order)))
            %% This makes the order flat
            %              test = xy_pair(center,1,1).*ones(1,size(order,2));
            %              xy_pair(center,:,2) = polyval(order_coeff,test);
            
            
            xy_pair(center,:,2) = polyval(order_coeff,xy_pair(center,:,1));
            
            % multiply by 100 to convert x,y mm to pixels then move 0,0 from center to bottom left
            xy_pair =xy_pair*100+4096/2;
            
            delta_w = diff(wavelength)/2;
            delta_w = [delta_w delta_w(end)];
            
            edges_w = [wavelength(1)-delta_w(1) wavelength+delta_w];
            if cheby ==1
                edges_x = interp1( finewavelengths, finex ,edges_w,'linear','extrap');
            else
                edges_x = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,edges_w,'linear','extrap');
            end
            %             edges_x = (-b + sqrt(b^2-4*a*(c-edges_w)))/(2*a);
            edges_x = edges_x*100+4096/2;
            
            
            %             offset = xy_pair(center,2:size(xy_pair,2),1)-xy_pair(center,1:size(xy_pair,2)-1,1);
            %             offset = [offset offset(end)];
            
            yoffset = [center-1:-1:1 0 -(1:size(order,1)-center)];
            
            detector = zeros(4096+4,4096+4); % 4096x4096 + 2 extra for padding
            
            %%% XY center pairs for entire order
            
            % xy_pair(vertical scalar index, horizontal vector of length, x and y)
            
            for jj = 1:size(xy_pair,2)
                for ii = [1:center-1 center+1:size(order,1)]
                    xy_pair(ii,jj,1) = xy_pair(center,jj,1);
                    %                     xy_pair(ii,jj,2) = xy_pair(center,jj,2)+ yoffset(ii)*offset(1);
                    xy_pair(ii,jj,2) = xy_pair(center,jj,2)+ yoffset(ii)*mean(diff(edges_x));
                    
                end
            end
            
            
            for ii = [1:center-1 center center+1:size(order,1)]
                
                for jj =1:size(xy_pair,2)
                    %                     xlist{ii}{jj} = xy_pair(ii,jj,1) + [-1 1 1 -1 -1]*offset(jj)/2 ; %x starting at top left, clockwise
                    %                     ylist{ii}{jj} = xy_pair(ii,jj,2) + [1 1 -1 -1  1]*offset(1)/2; % y ***extra final point to draw box
                    xlist{ii}{jj} = [edges_x(jj) edges_x(jj+1) edges_x(jj+1) edges_x(jj) edges_x(jj)];
                    ylist{ii}{jj} = xy_pair(ii,jj,2) + [1 1 -1 -1  1]*mean(diff(edges_x))/2; % y ***extra final point to draw box
                    
                    % +1 to the x and y list to account for a 1 cell padd
                    % to remove edge clipping errors.
                    [area,ind]=polyfillaa(2+xlist{ii}{jj},2+ylist{ii}{jj},4+4096, 4+4096);
                    ind = ind+1;
                    detector(ind) = detector(ind)+order(ii,jj).*area./sum(area);
                    clear ind area
                end
                
            end
            
            
            detector = detector(3:end-2,3:end-2); % remove padding
            detector = detector';
            
            
            
            % Plotting - be careful with this
            %             figure(100);
            %             hold on
            %
            %             for jj = 1:200%size(order,2)
            %
            %                 for ii = [1:center-1 center center+1:size(order,1)]
            %
            %
            % %                     plot(xy_pair(:,jj,1),xy_pair(:,jj,2),'.k')
            %
            %                     plot(xlist{ii}{jj},ylist{ii}{jj},'-','color',colors{rem(jj,length(colors)-1)+1})
            %
            %
            %
            %                 end
            %             end
            
            
        end
        
        function WriteFits(filename,image,headerinfo)
            
            import matlab.io.*
            DataType = 'single';
            % create fits file
            fptr = fits.createFile(filename);
            
            % Create Image extension
            fits.createImg(fptr,DataType,size(image));
            
            % Write image to file
            fits.writeImg(fptr,image);
            
            % Write additional keywords
            for ii = 1:size(headerinfo,1)
                fits.writeKey(fptr,headerinfo{ii,1},headerinfo{ii,2},headerinfo{ii,3});
                
            end
            % close file
            fits.closeFile(fptr);
            
        end
        
        function [trimWave,trimCounts] = trimToBand(wavelength,counts,bandPass)
            
            ind = wavelength<=bandPass(2) & wavelength>=bandPass(1);
            
            trimWave = wavelength(ind);
            trimCounts = counts(ind);
            
            
        end
            
            
               
    end
end

function [areas,ret] = polyfillaa(px,py,sx,sy)



% Clip grid to the enclosing box

left = max(floor(min(px)),0);
right = min(floor(max(px)),sx-1);
bottom = max(floor(min(py)),0);
top = min(floor(max(py)),sy-1);
nx=right-left+1;
ny=top-bottom+1;


if nx < 1 || ny < 1
    
%     disp('clipping missed grid')
    ret = 0;
    areas = 0;
    return
end

%npix is the maximum possible number of clipped polys
npix=nx*ny;

if npix <= 0
%     disp('clipping missed grid')
    ret = 0;
    areas = 0;
    return
end

ind = 1;

areas = [];%zeros(1,npix);
ret = [];%zeros(1,npix);
% 
% figure(1)
% plot([px,px(1)],[py,py(1)],'-k')

for p = 1
    px1 = px;
    py1 = py;
    
    
    for j = bottom(p):top(p)
        for i =left(p):right(p)
            px_out = px1;
            py_out = py1;
%             polyclip(i,j,px_out,py_out)
            Polygon=[px_out;py_out]';
            unitPixel = [i,i+1,i+1,i;j,j,j+1,j+1]';
            
            clippedPolygon = sutherlandHodgman(Polygon,unitPixel);
%             figure(1)
%             hold on
                 
            
            if isempty(clippedPolygon) ==1
                continue 
                
            end
%             plot([clippedPolygon(:,1)',clippedPolygon(1,1)],[clippedPolygon(:,2)',clippedPolygon(1,2)],'-o')
           
            
            ret(ind)=i+j*sx;
            
            
            px_out = clippedPolygon(:,1);
            py_out = clippedPolygon(:,2);
            

                       
            areas(ind) = abs(sum(px_out.* circshift(py_out,-1)-py_out.*circshift(px_out,-1))./2);
            ind = ind+1;
        end
    end
end


zero_ind = find(areas==0);
areas(zero_ind)=[];
ret(zero_ind)=[];

end
function clippedPolygon = sutherlandHodgman(subjectPolygon,clipPolygon)
%The inputs are a table of x-y pairs for the verticies of the subject
%polygon and boundary polygon. (x values in column 1 and y values in column
%2) The output is a table of x-y pairs for the clipped version of the 
%subject polygon. 
%% Helper Functions
 
    %computerIntersection() assumes the two lines intersect
    function intersection = computeIntersection(line1,line2)
 
        %this is an implementation of
        %http://en.wikipedia.org/wiki/Line-line_intersection
 
        intersection = zeros(1,2);
 
        detL1 = det(line1);
        detL2 = det(line2);
 
        detL1x = det([line1(:,1),[1;1]]);
        detL1y = det([line1(:,2),[1;1]]);
 
        detL2x = det([line2(:,1),[1;1]]);
        detL2y = det([line2(:,2),[1;1]]);
 
        denominator = det([detL1x detL1y;detL2x detL2y]);
 
        intersection(1) = det([detL1 detL1x;detL2 detL2x]) / denominator;
        intersection(2) = det([detL1 detL1y;detL2 detL2y]) / denominator;
 
    end %computeIntersection
 
    %inside() assumes the boundary is oriented counter-clockwise
    function in = inside(point,boundary)
 
        pointPositionVector = [diff([point;boundary(1,:)]) 0];
        boundaryVector = [diff(boundary) 0];
        crossVector = cross(pointPositionVector,boundaryVector);
 
        if ( crossVector(3) <= 0 )
            in = true;
        else
            in = false;
        end
 
    end %inside
 
%% Sutherland-Hodgman Algorithm
 
    clippedPolygon = subjectPolygon;
    numVerticies = size(clipPolygon,1);
    clipVertexPrevious = clipPolygon(end,:);
 
    for clipVertex = (1:numVerticies)
 
        clipBoundary = [clipPolygon(clipVertex,:) ; clipVertexPrevious];
 
        inputList = clippedPolygon;
 
        clippedPolygon = [];
        if ~isempty(inputList)
            previousVertex = inputList(end,:);
        end
 
        for subjectVertex = (1:size(inputList,1))
 
            if ( inside(inputList(subjectVertex,:),clipBoundary) )
 
                if( not(inside(previousVertex,clipBoundary)) )  
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);
                end
 
                clippedPolygon(end+1,1:2) = inputList(subjectVertex,:);
 
            elseif( inside(previousVertex,clipBoundary) )
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);                            
            end
 
            previousVertex = inputList(subjectVertex,:);
            clipVertexPrevious = clipPolygon(clipVertex,:);
 
        end %for subject verticies                
    end %for boundary verticies
end %sutherlandHodgman
function F=circ_gauss(X,Y,Sigma,center)
%--------------------------------------------------------------------------
% circ_gauss function                                                General
% Description: Calculate 2D circular Gaussian in a 2-D grid.
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D Gaussian.
%          - same as the x-ccordinates, but for the y-axis.
%          - Sigma of the Gaussian or [SigmaX, SigmaY] in case sigma
%            is different for each axis.
%            By default SigmaY=SigmaX.
%            If empty matrix use default.Si
%          - Center of the Gaussian [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - Maximum radius of Gaussian behond to set it to zero.
%            Default is Inf.
%            MaxRad is measured from the center of the kernel and not
%            the center of the Gaussian.


% Example: 
%          F=circ_gauss(MatX,MatY,[1],[0 0]);
%          surface(F);
%--------------------------------------------------------------------------



SigmaX = Sigma(1);
SigmaY = Sigma(2);

X0 = center(1);
Y0 = center(2);

F = 1./(2.*pi.*SigmaX.*SigmaY) .*exp(-1./(2.).* ((X-X0).^2./SigmaX.^2 +(Y-Y0).^2./SigmaY.^2));

F = F./sum(sum(F));

% set elements outside MaxRad to zero:
% if (~isinf(cutoff)),
%    MatR = sqrt(X.^2 + Y.^2);
%    I = find(MatR>cutoff);
%    F(I) = 0;
% end
% 
% if (isnan(Norm)),
%    % do not normalize
% else
%    F = Norm.*F./sumnd(F);
% end
end    
    
    
