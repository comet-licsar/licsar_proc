function [out, frfreqx, frfreqy] = cpxdetremd(data,numfringesx,numfringesy)
% CPXDETREND -- remove phase trend from complex data
%   Trend is notated in fringes over the image, not in frequency.
%   It is complex subtracted, for real valued phase DATA use DETREND,
%   or use DATA=complex(cos(DATA),sin(DATA)); (but this wrappes DATA).
%
%   CPXDETREND simulates complex interferogram with some noisy
%   fringes (DATA), estimates fringerates, subtract trend, and
%   plot result for demonstative purposes.
%  
%   CPXDETREND (DATA) returns detrended DATA. The number of fringes
%   over the image is estimated, shown, and this trend is removed.
%   Estimate is over total image by peak
%   estimation in Fourier domain (sum over all lines, columns).
%   Estimated is the number of integer fringes over the domain of,
%   the data; (no oversampling of data is performed). 
%   Fringe frequency is not given, simply the number of integer
%   fringes over the data in 2 drections.
%   DATA is assumed to be complex data, stored in major row order,
%   pixel interleaved (mph) format.
%
%   [OUT, FRX, FRY] = CPXDETREND (DATA) optionally outputs the 
%   Fourier transforms where the peak estimation is performed.
%
%   CPXDETREND (DATA, NUMFRINGESX, NUMFRINGESY) removes NUMFRINGESX
%   in X direction and NUMFRINGESY in Y direction.
%   A negative number of fringes indicates addition of an upward trend.
%
%   Large trends cannot be estimated with cpxdetrend.
%
%   Examples:
%     To remove 5.3 fringes in x direction, use:
%       OUT = CPXDETREND(DATA,5.3,0);
%       imagesc(angle(OUT));
%
%     To estimate the number of fringes, remove it, and check
%     the peak estimation, use:
%       [OUT, FRX, FRY] = CPXDETREND(DATA);
%     or, if you want to use the simulated data for testing:
%       [OUT, FRX, FRY] = cpxdetrend;
%       figure(1);
%       subplot(2,1,1), plot(FRX); title ('Peak in X direction');
%       subplot(2,1,2), plot(FRY); title ('Peak in Y direction');
%
%   See also: FREADBK, SIMULATESLC, FFT, DETREND, OVERSAMPLE
%

% $Revision: 1.5 $  $Date: 2001/09/28 14:24:31 $
% Bert Kampes, 27-Jun-2000

% probably better to zeropadd input so that spectrum is oversampled
% to estimate half fringes
%// BK 07-Aug-2001

%%% Check input
twopi    = 2. .* pi;
needhelp = 1;
simulate = 0;
estimate = 0;
if     (nargin==0)
  disp('Simulating data...');
  noiselevel = 1.2;
  needhelp = 0;
  simulate = 1;
  estimate = 1;
  L        = 100; 
  P        = 80;
  numfringesx = 5.5;%			simulate 5.5 fringes in x direction
  numfringesy = 3.8;%			simulate 3.8 fringes in y direction
  dx       = twopi ./ P;
  trendx   = 0:dx:twopi-dx;
  trendx   = (numfringesx) .* (trendx-pi);
  dy       = twopi ./ L;
  trendy   = 0:dy:twopi-dy;
  trendy   = (numfringesy) .* (trendy-pi);
  data     = ones(L,1) * lying(trendx) + standing(trendy) * ones(1,P);
  data     = complex(cos(data),sin(data));
  data     = complex(real(data)+noiselevel.*randn(L,P), ...
		     imag(data)+noiselevel.*randn(L,P));
  %trendxy  = (1./(L*P)) .* standing(trendy) * lying(trendx);
  %data     = complex(cos(trendxy),sin(trendxy));
elseif (nargin==1)
  needhelp = 0;
  estimate = 1;
elseif (nargin==3)
  disp('removing trend...');
  needhelp = 0;
end
if (isreal(data)) needhelp=1; end;
if (needhelp) helphelp; end;


%%% variables
[L P] = size(data);%		L lines (Y); P pixels (X)

%%% compute fringefreq. how to estimate half a fringe???
% fft over all lines seems like overkill...
if (estimate) 
  %warning('not ok yet, better oversample...');
  disp('estimating fringerates by FFT...');
  % x direction
  frfreqx = sum((abs(fft(data,[],2))),1);
  [maxval,ii] = max(frfreqx);
  SNRx = P .* maxval ./ sum(frfreqx);
  disp(['SNRx: ']);
  disp(SNRx);
  numfringesx = ii - 1;%		matlab starts array at 1
  if ( numfringesx > P/2 )
    numfringesx = numfringesx - P;%		negative number
  end;
  disp(['numfringesx: ']);
  disp(numfringesx);
  disp(' fringes over ');
  disp(P);
  disp(' pixels = ');
  disp(twopi*numfringesx/P);
  disp(' rad/pixel');

  % y direction
  frfreqy = sum((abs(fft(data,[],1))),2);
  [maxval,ii] = max(frfreqy);
  SNRy = L .* maxval ./ sum(frfreqy);
  disp(['SNRy: ']);
  disp(SNRy);
  numfringesy = ii - 1;%		matlab starts array at 1
  if ( numfringesy > L/2 )
    numfringesy = numfringesy - L;%		negative number
  end;
  disp(['numfringesy: ']);
  disp(numfringesy);
  disp(' fringes over ');
  disp(L);
  disp(' pixels = ');
  disp(twopi*numfringesx/L);
  disp(' rad/pixel');
end

%%% do the subtraction complex notation
disp(['subtracting fringes: ']);
disp(numfringesx);
disp(numfringesy);
dx     = twopi / P;
trendx = 0:dx:twopi-dx;
trendx = lying(numfringesx .* trendx);
dy     = twopi / L;
trendy = 0:dy:twopi-dy;
trendy = standing(numfringesy .* trendy);

disp('here1');
% remove trend line by line for memory considerations
trendx = complex(cos(trendx),-sin(trendx));
trendy = complex(cos(trendy),-sin(trendy));
out    = zeros(L,P);
for ii=1:L
  out(ii,:) = data(ii,:) .* trendx;
end
disp('here2');
for ii=1:P
  out(:,ii) = out(:,ii) .* trendy;
end
disp('here3');

