function [fig_handle,h1,h2] = E_DTI_imagescn_overlay_5(I, I2, mask_I, mask_I2, scale, VDims, dims, afd_base,FigureWidth, timedimension)
% function imagescn(I,[min max],[rows cols],FigureWidth,timedimension)
%
% function to display multiple images
%	I can be 2-d, i.e., single image
%            3-d,       array of images
%            4-d,       2-d array of images
%            5-d,       3-d array of images
%
%   timedimension specifies which dimension is used for time in "movie" mode
%        and activates the movie mode button. if time dimension is not
%        specified then there are 5 dimensions and the last dimension
%        defaults to the timefenghua
%        dimension. the time dimension must be between 3 and number of
%        dimensions (max=5).
%
% user specified scale [min max] or is applied to all images (i.e., clims)
% or scale=[min1 max1; min2 max2; ...] will be applied to images 1,2,...,etc.
% defaults to independent scaling of each sub-image to full range (for scale=[])
%
% sub-images are displayed in N1 rows by N2 columns
%   N1 and N2 may be specified in function input as "rows" and "cols"
%		if input values for rows & cols input are such that
%		rows*cols < full number of images, then the 1-st rows*cols images are displayed
%	or N1 and N2 may be calculated by default:
%		for 3-d array (N1 and N2 are calculated to make approx square image array)
%		for 4-d array (N1 and N2 are 3rd and 4th array dimensions, respectively)
%
% FigureWidth sets the width in inches (defaults to 6 inches). It also sets the paper width
% so that exported figures (to files) will have the specified width.
%
% usage:  imagescn(I)
%         imagescn(I,[],[],[],[])
%         imagescn(I,scale)
%         imagescn(I,[],[rows cols])
%         imagescn(I,scale,[rows cols])
%         imagescn(I,[],[],[],timedimension)
%         ...

% written by: 	Peter Kellman  (kellman@nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%   tools by: 	Dan Herzka  (herzkad@nhlbi.nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%
Nd=ndims(I);
if ~exist('timedimension')& Nd==5; timedimension=5; moviemodeflag=1; end 
if ~exist('timedimension'); timedimension=[]; moviemodeflag=0; end
if nargin==5
    if isempty(timedimension); moviemodeflag=0; else; moviemodeflag=1;end
    if timedimension>Nd | timedimension<=2; moviemodeflag=0;end
    if moviemodeflag==1
        if Nd==4 & timedimension==3
            I=permute(I,[1 2 4 3]);
        end
        if Nd==5 & timedimension==3
            I=permute(I,[1 2 4 5 3]);
        end
        if Nd==5 & timedimension==4
            I=permute(I,[1 2 3 5 4]);
        end
    end
end
ntimeframes=size(I,Nd);
if Nd==2; % case of single image
    N=1;
	N1=1; N2=1;
elseif Nd==3; % case of array of images
    if moviemodeflag==0
        N=size(I,3);
		N2=ceil(sqrt(N)); N1=ceil(N/N2);
	else
        N=1;N1=1;N2=1;
	end
elseif Nd==4; % case of 2-d array of images
    if moviemodeflag==0
        N1=size(I,3);
        N2=size(I,4);
		N=N1*N2;
	else
        N=size(I,3);
		N2=ceil(sqrt(N)); N1=ceil(N/N2);
	end
elseif moviemodeflag==1 & Nd==5; % case of 2-d array of images
    N1=size(I,3);
    N2=size(I,4);
	N=N1*N2;
end
if exist('dims','var');
	if length(dims)==2; rows=dims(1);cols=dims(2);
		N1=rows;N2=cols;
	else
		if ~isempty(dims);disp('Error: must enter [rows cols] for dimensions'); return;end
	end
end
if ~exist('scale','var'); scale=[];end
if 1<size(scale,1) & size(scale,1) < min(N,N1*N2)
    disp('scale vector must be either: empty, size 1x2, i.e. [min max],')
    disp('or a matrix [min1 max1;min2 max2;...] with # rows equal to number of plots');
    return
end
                
set(0,'Units','Inches');
scnsize=get(0,'ScreenSize'); % [left,bottom,width,height] % full screen size available
ScreenWidth=scnsize(3); ScreenHeight=scnsize(4); % width & height in inches
Xsize=size(I,2);Ysize=size(I,1); % size of pixels in image (Xsize x Ysize)
border_percent=.03;
deltaX =(border_percent*Xsize); deltaY=(border_percent*Ysize); % calculate pixel borders as % of size
X0size=N2*Xsize+(N2-1)*deltaX; Y0size=N1*Ysize+(N1-1)*deltaY; % full figure size in pixels (before display)
aspect_ratio=Y0size/X0size; % aspect ratio

% center figure on screen with specified figure width in inches
if ~exist('FigureWidth'); FigureWidth=6;end  % figure width in inches (default is 6")
if isempty(FigureWidth); FigureWidth=6;end
FigureHeight=FigureWidth*aspect_ratio; % figure height in inches
FigureBottom=(ScreenHeight-FigureHeight)/2;
FigureLeft=(ScreenWidth-FigureWidth)/2;
fig_handle=figure('NumberTitle','off','color','k','Menubar','figure','Visible','off');
E_DTI_Set_Fig_Icon(fig_handle);
set(fig_handle,'Visible','on')
set(fig_handle,'Units','Inches')
set(fig_handle,'Position',[FigureLeft FigureBottom FigureWidth FigureHeight])

% calculate sub-image dimensions in inches
SubImageWidth=FigureWidth*Xsize/X0size;
SubImageHeight=FigureHeight*Ysize/Y0size;
Xborder=FigureWidth*deltaX/X0size;
Yborder=FigureHeight*deltaY/Y0size;

% set background color to be white
set(fig_handle,'Color',[0 0 0]);

% calculate sub-image dimensions in normalized units
SubImageWidth=SubImageWidth/FigureWidth;
SubImageHeight=SubImageHeight/FigureHeight;
Xborder=Xborder/FigureWidth;
Yborder=Yborder/FigureHeight;



cmin = afd_base(1,:); % max(0,min(I2(:)));
cmax = afd_base(2,:); % min(100,max(I2(:)));
mask_I2 = mask_I2>0;

m = 64;
m2 = 64; %ceil(cmax)-floor(cmin)+1;

cmap2 = {parula(m2),pink(m2),pink(m2),parula(m2),parula(m2),pink(m2),pink(m2),parula(m2),parula(m2),parula(m2)};

% cmap = [gray(64);parula(64)];
%
for k=1:min(N,N1*N2)
	i=ceil(k/N2);
	j=k-(i-1)*N2;
    if Nd>3
    	i0=ceil(k/size(I,4));
		j0=k-(i0-1)*size(I,4);
    end        
	SubImageLeft=(j-1)*(SubImageWidth+Xborder);
	SubImageBottom=(N1-i)*(SubImageHeight+Yborder);
	subplot('Position',[SubImageLeft SubImageBottom SubImageWidth SubImageHeight])
	if isempty(scale)
 %       
		h1 = imagesc(I(:,:,k)); axis image; axis off; hold on
        c1 = get(h1,'CData');
        
        cmin1 = min(c1(:));
        cmax1 = max(c1(:));
        C1 = min(m,round((m-1)*(c1-cmin1)/(cmax1-cmin1))+1);
        C1(~mask_I(:,:,k)) = 1;
        
        h2 = imagesc(I2(:,:,k),'AlphaData',mask_I2(:,:,k)); axis image; axis off; hold off
        c2 = get(h2,'CData'); 
        

        cmin2 = max(cmin(k)); %,round(min(c2(mask_I2(:,:,k))))); %min(c2(:))
        cmax2 = min(cmax(k)); %,round(max(c2(mask_I2(:,:,k))))); %max(c2(:))
        C2 = min(m2,round((m2-1)*(c2-cmin2)/(cmax2-cmin2))+1);
        C2(mask_I2(:,:,k)) = m+C2(mask_I2(:,:,k));
        C2(~mask_I2(:,:,k)) = C1(~mask_I2(:,:,k));
        C2(c2<cmin(k)) = C1(c2<cmin(k)); 
%         C2(c2<cmin(k)) = m+m2+1; %C2(c2>32) = m+m2+1;
        disp(['The median is ' num2str(median(c2(mask_I2(:,:,k)))) ' ' num2str(median(C2(mask_I2(:,:,k))))])
%         disp(['The median is ' ])
        
        set(h1,'CData',C1);
        set(h2,'CData',C2,'AlphaData',mask_I2(:,:,k));
        
        cmap = [gray(m);cmap2{1,k};cool(1)]; %[gray(m);parula(m2)];
%         caxis([min(C1(:)) max(C2(:))])
        caxis([min(C1(:)) m+m2+1])
        colormap(gca,cmap)
%         if k==1
%             disp(['cmin2 and cmax2 is [' num2str(round(min(c2(mask_I2(:,:,k))))) ' ' num2str(round(max(c2(mask_I2(:,:,k))))) ']'])
%             disp(['cmin2 and cmax2 is [' num2str(cmin2) ' ' num2str(cmax2) ']'])
%         end
%         h3 = colorbar;
%         set(h3,'YLim',[m+1 m+m2])
%         set(h3,'Ticks',median(C2(mask_I2(:,:,k))))
  %      
    elseif size(scale,1)==1
		h1 = imagesc(I(:,:,k),scale); axis image;axis off; hold on
        c1 = get(h1,'CData');
        
        cmin1 = min(c1(:));
        cmax1 = max(c1(:));
        C1 = min(m,round((m-1)*(c1-cmin1)/(cmax1-cmin1))+1);
        C1(~mask_I(:,:,k)) = 1;

        h2 = imagesc(I2(:,:,k),'AlphaData',mask_I2(:,:,k)); axis image; axis off; hold off
        c2 = get(h2,'CData');

        cmin2 = cmin(k); %min(c2(:));
        cmax2 = cmax(k); %max(c2(:));
        C2 = min(m2,round((m2-1)*(c2-cmin2)/(cmax2-cmin2))+1);
        C2(mask_I2(:,:,k)) = m+C2(mask_I2(:,:,k));
        C2(~mask_I2(:,:,k)) = C1(~mask_I2(:,:,k));
        C2(c2<cmin(k)) = m+1; %C1(c2<cmin(k));  if shown as c2(c2<cmin) = cmin;

        set(h2,'CData',C2,'AlphaData',mask_I2(:,:,k));
        set(h1,'CData',C1);
        
        cmap = [gray(m);parula(m2)];
%         caxis([min(C1(:)) max(C2(:))])
        caxis([min(C1(:)) m+m2])
        colormap(cmap)
        
%         if k==1
%             disp(['cmin and cmax is [' num2str(cmin) ' ' num2str(cmax) ']'])
%         end
        
    elseif size(scale,1)==min(N,N1*N2);
        imagesc(I(:,:,k),scale(k,:)); axis image;axis off
    end
    if exist('VDims')
    set(gca,'DataAspectRatio',VDims)
    end
    if moviemodeflag==1        
        if Nd==3
        	setappdata(gca, 'ImageData', I);
        elseif Nd==4
        	setappdata(gca, 'ImageData', squeeze(I(:,:,k,:)));
        elseif Nd==5
            setappdata(gca, 'ImageData', squeeze(I(:,:,i0,j0,:)));
        end
		setappdata(gca, 'ImageRange', [1 ntimeframes]);
		setappdata(gca, 'ImageRangeAll', [1 ntimeframes]);
		setappdata(gca, 'CurrentImage', 1);
    end
end

set(fig_handle, 'PaperPosition', [0 0 FigureWidth FigureHeight]);


% WL_tool  % launch window level tool (button)
% PZ_tool  % launch Pan-Zoom tool (button)
% ROI_tool % launch ROI tool (button)
% MV_tool  % launch movie tool (button)
% PM_tool  % launch pixel (point) measurement tool (button)
% RT_tool  % launch rotate/flip tool (button)

% colormap(gray)
% colormap(parula)
% colorbar

