function [avg,nstd] = robmean(x,throw)
% function [avg,nstd] = robmean(x,throw)
% "robust" mean: removes NaNs and removes highest and lowest values (top
% throw % AND bottom throw %) before calculating mean


if min(size(x))==1
  bad = find(~isfinite(x));
  x(bad) = [];
  
  drop = round(throw*length(x));
  sx = sort(x);
  avg = mean(sx(drop+1:length(x)-drop));
  if nargout==2
    nstd = std(sx(drop+1:length(x)-drop));
  end;
else
  avg = nan*ones(1,size(x,2));
  if nargout==2
    nstd = avg;
  end;
  for i = 1:size(x,2),
    xc = x(:,i);
    bad = find(~isfinite(xc));
    xc(bad) = [];
    
    drop = round(throw*length(xc));
    sxc = sort(xc);
    avg(i) = mean(sxc(drop+1:length(xc)-drop));
  if nargout==2
    nstd(i) = std(sxc(drop+1:length(xc)-drop));
  end;
  end; % for i = 1:size(x,2),
end % if min(size(x))==1
