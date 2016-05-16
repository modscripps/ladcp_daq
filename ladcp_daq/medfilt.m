function [xp] = medfilt(x,wlen);
% MEDFILT  Median filter.
%   Applies a running window to a data series and replaces each point with
%   the median of the values in window centered on that point.
%

if nargin==1
  wlen = 3;
end
if mod(wlen,2)==0
  % wlen must be ODD
  wlen = wlen-1;
end

[m,n] = size(x);
if m>n
  % x needs to be a row vector
  colvec = 1;
  x = x';
  n = m;
else
  colvec = 0;
end
halfwin = (wlen-1)/2;
xp = x;

for i = 1:n,
  if i==1
    % Beginning
    %xw = [x(i:i+halfwin)];
    xw = x(i);
  elseif i==n
    % End
    %xw = [x(i-halfwin:i)];
    xw = x(i);
  elseif i<=halfwin
    % Beginning segment
    xw = [x(1:i+halfwin)];
  elseif i>n-halfwin
    % End segment
    xw = [x(i-halfwin:n)];
  else
    % Middle segment
    xw = [x(i-halfwin:i+halfwin)];
  end
  xwm = nanmedian(xw);
  % Replace with mean
  xp(i) = xwm;
end
  
if colvec
  xp = xp';
end

%ind = 1:n;
%plot(ind,x,'*',ind,xp,'r');

return
