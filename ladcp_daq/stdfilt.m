function [xp] = stdfilt(x,wlen,nsig);
% STDFILT  Standard deviation filter.
%   Applies a running window to a data series and removes outliers based on a
%   standard deviation criterion (outside of +/- nsig std devs)
%
% A few choices: - Could remove a linear or quadratic fit from each window first
%   - Could replace outliers with mean value within window or interpolate
%   - Probably should just consider center point of window and don't include
%      its value in the std calculation
%   - At ends, need to use a lopsided window

if nargin==1
  wlen = 20;
end
if nargin<3
  nsig = 2;
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
    xw = [x(i+1):x(i+halfwin)];
  elseif i==n
    % End
    xw = [x(i-halfwin:i-1)];
  elseif i<=halfwin
    % Beginning segment
    xw = [x(1:i-1) x(i+1:i+halfwin)];
  elseif i>n-halfwin
    % End segment
    xw = [x(i-halfwin:i-1) x(i+1:n)];
  else
    % Middle segment
    xw = [x(i-halfwin:i-1) x(i+1:i+halfwin)];
  end
  xwm = nanmean(xw);
  xwsig = nanstd(xw);
  % Go straight to std of x over window
  if abs(x(i)-xwm) > nsig*xwsig
    % Replace with mean
    xp(i) = xwm;
  end
end
  
if colvec
  xp = xp';
end

%ind = 1:n;
%plot(ind,x,'*',ind,xp,'r');

return
