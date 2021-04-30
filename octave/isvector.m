function res = isvector(in)
% ISVECTOR  --  True for vectors, ie, if "min(size(a))==1".
%
%   Only for 2 dimensional matrices/vectors/scalars.
%   If IN is a scalar then it is not a vector...
%

%// BK 07-Aug-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:45 $

res=0;
if (length(size(in))~=2)
  warning('not intented for these things.');
elseif (min(size(in))==1 & ~isscalar(in))
  res=1;
end

%%% EOF
