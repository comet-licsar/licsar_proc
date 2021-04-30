function res = isscalar(in)
% ISSCALAR  --  True for scalars, ie, if "max(size(a))==1".
%
%   Mainly used to test input in functions.
%

%// BK 07-Aug-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:45 $

res=(max(size(in))==1);

%%% EOF
