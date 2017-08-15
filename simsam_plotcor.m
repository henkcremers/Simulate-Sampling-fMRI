function [scatfig,s] = simsam_plotcor(data,varargin)
filename = 'corrplot';
saveplot = 0; 
fsize = 35; % font size
psize = 10; 

% Input.
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'saveplot'
                saveplot = 1;
                filename = varargin{i+1};
            case 'fsize'
                fsize = varargin{i+1};
            case 'psize'
                psize = varargin{i+1};    
        end
    end
end

V1 = data(:,1); 
V2 = data(:,2);
c = corr(V1,V2);

% set the location depending on correlation value
if c>0; yloc = 0.82; elseif c<0; yloc = 0.18; end


%% create the figure
scatfig = figure;
legtxt = ['r = ' num2str(double(c))]; legtxt = legtxt(1:9);
s = scatter(V1,V2);
set(s,'MarkerFaceColor','b');

% add the correlation coef
annotation('textbox',...
    [0.15 yloc 0.3 0.1],...
    'String',legtxt,...
    'FontSize',fsize,...
    'FontName','Arial',...
    'LineWidth',2,...
    'BackgroundColor','w'); 

%% plot a regression line + confidence interval

[coef_fit,s] = polyfit(V1,V2,1);
hold on
[yfit,dy] = polyconf(coef_fit,xlim,s,'predopt','curve');
plot(xlim,yfit,'color','k','LineWidth',4)
line(xlim,yfit-dy,'color','k','linestyle',':','LineWidth',2)
line(xlim,yfit+dy,'color','k','linestyle',':','LineWidth',2)
xl = xlabel('Behavioral Variable','FontSize',fsize);
yl = ylabel('Brain Activation','FontSize',fsize);

% fix the background
set(gca,'Color',[0.9 0.9 0.9]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
set(0,'DefaultFigureColor',[1 1 1])
%set(get(gca,'xlabel'),'position',[-0.05 0 0])
set(gca,'Fontname','Gill Sans')

if saveplot == 1
axis tight
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'PaperPosition', [0 0 psize psize]);
set(gcf, 'PaperSize', [psize psize]);         
filename = [filename];
saveas(gcf,filename,'pdf');
%saveas(gcf,filename,'bmp');
end

return