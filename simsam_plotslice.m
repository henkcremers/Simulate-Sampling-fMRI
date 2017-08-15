function [braindata,blobcontour,underlay] = simsam_plotslice(braindata,ctresh,varargin)

% defaults
orientation = 'hor';
slicenum = 35;
underlay = [];
plotcir = 0;
contour = [];
saveplot = 0;
smoothdat = 0;
filename = 'brainplot';

% get the input.
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'or'
                orientation = varargin{i+1};
            case 'slicenum'
                slicenum = varargin{i+1};
            case 'underlay'
                underlay = varargin{i+1};
            case 'contour'
                contour = varargin{i+1};
            case 'plotcir'
                plotcir = 1;
                voxplot = varargin{i+1};
            case 'smoothdat'
                smoothdat = 1;
                %voxplot = varargin{i+1};    
            case 'saveplot'
                saveplot = 1;
                filename = varargin{i+1};
            case 'name'    
                filename = varargin{i+1};
        end
    end
end


if isempty(underlay)
    
    % only needed when there is no bgdata - load the SPM grey matter data
    load spmGMdat.mat
    Data  = GM.Data;
    coord = GM.coord;
    
    % set the coordinates and frame.
    x = unique(coord(:,1));
    y = unique(coord(:,2));
    z = unique(coord(:,3));
    
    nx = length(x);
    ny = length(y);
    nz = length(z);
    
    switch orientation
        
        case 'hor'
            % % slice coord - use horizontal plane
            zslice = coord(:,3)==slicenum;
            slicecoord = coord(zslice,1:2);
            voxelval = Data(:,zslice);
            slicedat = accumarray([slicecoord(:,1),slicecoord(:,2)],voxelval,[max(x) max(y)]);
            yaxis = 1:max(y);
            xaxis = 1:max(x);
            % crop the image
            braindata = braindata(1:max(y),1:max(x));
        case 'sag'
            xslice = coord(:,1)==slicenum;
            slicecoord = coord(xslice,2:3);
            voxelval = Data(:,xslice);
            slicedat = accumarray([slicecoord(:,1),slicecoord(:,2)],voxelval,[max(y) max(z)]);
            yaxis = 1:max(z);
            xaxis = 1:max(y);
            % crop the image
            braindata = braindata(1:max(z),1:max(y));
        case 'cor'
            yslice = coord(:,2)==slicenum;
            slicecoord = coord(yslice,[1 3]);
            voxelval = Data(:,yslice);
            slicedat = accumarray([slicecoord(:,1),slicecoord(:,2)],voxelval,[max(x) max(z)]);
            yaxis = 1:max(z);
            xaxis = 1:max(x);
            % crop the image
            braindata = braindata(1:max(z),1:max(x));
    end

    % fix the rotation of the background
    bgdat = rot90(slicedat,3);
    
    % create the brain contour data.
    braincontour = zeros(size(bgdat,1),size(bgdat,2));
    inbrain = (bgdat>0.05);
    braincontour(inbrain) = 1;
    
    % use 2 colormaps - underlay background.
    mx = max(max(bgdat));
    mi = min(min(bgdat));
    rg = mx-mi;
    sfac = rg/2;
    bgdat = bgdat./sfac;
    tdown = max(max(bgdat)) + 1;
    bgdat = bgdat -tdown;
    
    % sumarize in structure
    underlay.bgdat = bgdat;
    underlay.braincontour = braincontour;
    underlay.inbrain = inbrain;
    
else
    bgdat = underlay.bgdat;
    braincontour= underlay.braincontour;
    inbrain =  underlay.inbrain;
    xaxis = 1:size(bgdat,2);
    yaxis = 1:size(bgdat,1);
end
%% 

% blob contour
if isempty(contour);
    blobcontour = (braindata<=-ctresh | braindata>=ctresh); 
else
    blobcontour = contour;
end

% use bg colors.
bgloc = braindata>-ctresh & braindata<ctresh;
brainplotdata = braindata;
brainplotdata(bgloc) = bgdat(bgloc);
% white background - set to 1
brainplotdata(~inbrain) = 1; 

%% create the figure
figure('Name',filename,'Position', [50, 100, 480, 600])
hold on

% contour of blobs - white
[dnj bc] = imcontour(xaxis,yaxis,blobcontour);
set(bc,'LineColor','w','LineWidth',0.1)

% plot the brain and underlay data
imagesc(xaxis,yaxis,brainplotdata);

% use the "double colormap"
load ('simsam_colmap_comb.mat');
caxis([-3,1])
colormap(combcmap)

% % contour of blobs - white
[dnj bc] = imcontour(xaxis,yaxis,blobcontour);
set(bc,'LineColor','w','LineWidth',0.1)

% contour of the brain - black
[dnj h] = imcontour(xaxis,yaxis,braincontour);
set(h,'LineColor','k','LineWidth',0.1)

if plotcir == 1;
    v1 = voxplot(1);
    v2 = voxplot(2);
    r = 2;
    ang=0:0.001:2*pi;
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(v1+xp,v2+yp,'g','LineWidth',3);
end

% remove the axis
set(gca,'visible','off')

% save the plot (TIFF/BMP)
if saveplot == 1
    axis tight
    imname = [filename  '.BMP']; %BMP  
    f=getframe(gca);
    [X, map] = frame2im(f);
    I2 = imcrop(X,[8, 38, 330, 420]);
    imwrite(I2,[imname],'BMP')   
end

return
