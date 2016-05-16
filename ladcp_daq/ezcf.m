function ax=ezcf(x,y,z,clim,ncolors);
%function ax=ezcf(x,y,z,clim,ncolors);
%4/09 MHA
%Make a pcolor-like plot with contourf.
if nargin < 4
    ncolors=32;
end
%'peg' values out of range
z(find(z<clim(1)))=clim(1);
z(find(z>clim(2)))=clim(2);

cvals=linspace(clim(1),clim(2),ncolors);
[c1,h]=contourf(x,y,z,cvals);
set(h,'edgecolor','none');
axis ij
caxis(clim)