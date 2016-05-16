function [zg,stdev,sterr] = gauss_sm1d(xd,zd,xg,scale)
% Smooth irregularly-spaced 1-d data with a Gaussian filter
  
bad = find(isnan(xd)|isnan(zd));
xd(bad) = [];
zd(bad) = [];

ndata = length(xd) ;
ngrid = length(xg) ;
zg = xg*NaN;
stdev = xg*NaN;
sterr = xg*NaN;

for k=1:ngrid
	dist   =  (xd-xg(k)).^2;
	weight =  exp( -dist / (scale*scale) ) ;
	wmax = max(weight);
	zg(k)  =  sum(zd.*weight) / sum(weight) ;
	stdev(k)  =  sqrt(sum((zd-zg(k)).^2.*weight) / sum(weight)) ;
	sterr(k)  =  sqrt(wmax*sum((zd-zg(k)).^2.*weight)) / sum(weight) ;
end

% Note: this routine uses my invented versions of standard deviation and
% standard error for a weighted average. Whether these quantities are
% meaningful (or could be defined better) is a matter of some debate...

return
