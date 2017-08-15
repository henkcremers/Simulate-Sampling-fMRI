%% Create an avi video of sub-sampling. 
% ====================================
close all; clear all; clc

%% settings
%  ======================================
ntot = 10000;  % full sample size
cval = 0.3;    % correlation values
std1 = 10;     % std variable 1 
std2 = 10;     % std variable 2 
m1 = 10;
m2 = 10;
mu = [m1 m2];
Sigma = [std1^2,cval*std1*std2;cval*std1*std2,std2^2];
fullsample = mvnrnd(mu,Sigma,ntot);
v2 = fullsample(:,2);
v1 = fullsample(:,1);
samples = 10:10:20;
ns = 10;
filename = 'subsampling_r=0.3.avi’;

% Prep video file.
% -----------------
vid = VideoWriter(filename);
vidObj.Quality = 200; 
vid.FrameRate = 7; 
bgcolor = [108/256 108/256 108/256];
or = [255/256 165/256 0];
open(vid);

csvals = []; 
c = [];
count = 0; 
pscat = [];
set(gca,'nextplot','replacechildren');

xlim = [min(v1) max(v1)];
ylim = [min(v2) max(v2)];

sigcount = 0;
totcount = 0; 
count = 0;

% Start
% ------

for n = samples;
    for k = 1:ns
    count = count + 1; 
    
    scatfig = figure;
    subplot(2,1,1);

    % layer 1 - full sample
    % --------------------
    s1 = scatter(v1,v2);
    
    if count==1; c = 'b'; else c = bgcolor; end
    
    set(s1,'MarkerFaceColor',c)
    set(s1,'MarkerEdgeColor',c)
    xlabel('Variable 1')
    ylabel('Variable 2')
    axis([min(v1) max(v1) min(v2) max(v2)])
    
    %plot a regression line..
    coef_fit1 = polyfit(v1,v2,1);
    y_fit1 = polyval(coef_fit1,xlim);
    hold on
    plot(xlim,y_fit1,'k');

    xlim = [min(v1) max(v1)];
    ylim = [min(v2) max(v2)];    
        
    % layer 2 - samples
    % --------------------
    population = 1:ntot;
    rsamp = randsample(population,n); 
    v2sample = v2(rsamp,1); 
    v1sample = v1(rsamp,1);

    [c pval] = corrcoef(v2sample,v1sample); c = c(2,1); pval = pval(2,1);

    if count>1
    s2 = scatter(v1sample,v2sample);
    set(s2,'MarkerFaceColor','b')
    set(s2,'MarkerEdgeColor','b')
    xlabel('Variable 1')
    ylabel('Variable 2')
    
    % plot a regression line..
    coef_fit2 = polyfit(v1sample,v2sample,1);
    y_fit2 = polyval(coef_fit2,xlim);
    hold on
    linecol = 'b';
    if pval<0.05; linecol = or; end
    p = plot(xlim,y_fit2);
    set(p,'Color',linecol);
    end
    
    subplot(2,1,2); 
    plot([min(samples) max(samples)],[0.3 0.3],'k') 
    xlabel('samples Size')
    ylabel('Correlation Coefficient')
    axis([min(samples) max(samples) -1 1])
    %axis([min(v1) max(v1) min(v2) max(v2)])
    hold on 
    
    xs = n;
    if count > 2
        csvals(count-1) = c;
        xsvals(count-1) = xs;
        
        %fix same correlation values...
        low = n-2;
        up = n+2;
        checkc = csvals(xsvals>low & xsvals<up);
        checkn = xsvals(xsvals>low & xsvals<up);
        samevals = sum(round(checkc,2)==round(c,2));
        if samevals > 0;
            sval = sign(randn(1));
            if sval >0
                maxx = max(checkn);
                xs = maxx + 0.1;
                if xs>up; xs = up-0.0001; end      
            elseif sval < 0
                minx = min(checkn);
                xs = minx - 0.1;
                if xs<low; xs = low+0.0001; end   
            end
        end
        
        %update values.. 
        csvals(count-1) = c;
        xsvals(count-1) = xs;
            
        scatter(xsvals,csvals,'b');
        axis([min(samples)-5 max(samples)+5 -1 1])
        %axis([min(v1) max(v1) min(v2) max(v2)])
        if pval < 0.05;
            sigcount = sigcount + 1;
            pscat1(sigcount) = xs;
            pscat2(sigcount) = c;
        end
        if sigcount >0         
            s3 = scatter(pscat1,pscat2,'filled');
            set(s3,'MarkerFaceColor',or)
            set(s3,'MarkerEdgeColor',or)
            axis([min(samples)-5 max(samples)+5 -1 1])
        end
    end

    % fix the background
    set(gca,'Color',[0.9 0.9 0.9]);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    set(0,'DefaultFigureColor',[1 1 1])
      
    % Write each frame to the file.
    F = getframe(scatfig);
    writeVideo(vid,F);
    close(scatfig)
    end
end

% Close the file.
close(vid)