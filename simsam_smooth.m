function [map_smooth] = simsam_smooth(map,sigma)
% repalce nan with 0
nanloc = isnan(map);
map(nanloc)=0;

% smooth the blobs
N = sigma.^2;
[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
map_smooth = conv2(map, f,'same'); 

% return the nan
map_smooth(nanloc)=NaN;
end

