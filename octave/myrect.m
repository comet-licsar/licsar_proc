function w = myrect(x)
% MYRECT  Weighting function for SAR spectra (range).
%    MYRECT(FR) computes rectangular window for FR (normalized).
%
%    This window is defined by:
%       rect(x) = 1,  abs(x) <= 0.5
%                 0,  otherwise
%
%    Note: rect not periodic (rect(x)==0 for all x>0.5).
%    To make it periodic with period [-1,1] simply change the input:
%    rect(x) -> rect(mod(x+1,2)-1))
%    This rescales x to the new interval: (-1,1], e.g.
%       x=-10:0.01:4;
%       plot(x,myrect(mod(x+1,2)-1));
%    To shift rect(x), simply change do rect(x-dx);
%    To scale rect(x), simply change do y=2*rect(x-dx);
%
%    This function is used in adaptive range filtering.
%    See also MYHAMMING, BOXCAR, HAMMING, RANGEDEMO, FILTRANGE, FILTERBLOCK.
%

% $Revision: 1.3 $  $Date: 2000/12/07 16:29:13 $
% Bert Kampes, 17/03/00

% Handle input
s = size(x);
if (max(s) ~= prod(s)) error('MYRECT: only scalars/vectors'); end;

% Compute rect for data
x=x(:);
w=zeros(size(x));
w(find(abs(x)<=.5))=1;
w=reshape(w,s);

%old
%w=zeros(size(x));
%w(find(abs(x)<=.5))=1;

%%% EOF
