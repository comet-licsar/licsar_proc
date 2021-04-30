function res = issamesize(a,varargin)
% ISSAMESIZE  --  True for arrays of the same size.
%
%   ISSAMESIZE(A,B) is true if all(size(a)==size(b))
%
%   ISSAMESIZE(A,B,C,...) returns 1 if all arrays have the same size.
%

%// BK 07-Aug-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:44 $

if (nargin<2) helphelp; break; end;
%%% Compare all input arrays with the first one.
for ii=1:length(varargin)
  res = (all(size(a)==size(varargin{ii})));
  if (res==0) return; end;
end
%res = all(size(a)==size(b));

%%% EOF
