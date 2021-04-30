function [n,scale,offset] = normalize(data,a,b)
%   NORMALIZE(DATA) normalizes data to interval [0,1].
% 
%   NORMALIZE(DATA,a,b) normalizes real DATA (in interval [c,d])
%   to interval [a,b].  Simple formula used:
%     DATA --> a + ((b-a)./(d-c)) .* (DATA-c);
%
%   [N, S, O] = NORMALIZE(...) optionally returns scale and offset
%   factors to normalize other values to the same grid.  Suppose DATA has
%   been normalized to interval [a,b], and now you want to transform X to 
%   the same interval:
%     [Y1, S,O] = NORMALIZE(DATA,a,b);
%      Y2       = S.*X + O;
%
%   See also INSAR toolbox.
%

%// $Revision: 1.3 $  $Date: 2001/09/28 14:24:32 $
%// Bert Kampes, 07-Feb-2001


%%% Set defaults.
if (~isreal(data))          error('Data should be real.'); end;
if (nargin~=1 & nargin~=3)  error('Not correct number input arguments'); end;
c = min(data(:));
d = max(data(:));

if (nargin == 1)
  %b = 1;
  %a = 0;
  scale  = 1./(d-c);
  offset = -c.*scale;
else
  offset = (d*a-c*b)./(d-c);
  scale  = (b-a)./(d-c);
end

%%% Normalize.
n      = scale.*data + offset;

%%% EOF.


