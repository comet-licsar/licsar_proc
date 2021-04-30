function dx = pyramid(nrows,ncolumns)
% PYRAMID --  Generate unit pyramid.
%   c = PYRAMID;        pyramid [0,1] of dimensions (512,512);
%   c = PYRAMID(L);     pyramid [0,1] of dimensions (L,L);
%   c = PYRAMID([L,P]); pyramid [0,1] of dimensions (L,P);
%   c = PYRAMID(L,P);   pyramid [0,1] of dimensions (L,P);
% 
%   PYRAMID is generated using find.
%
%   See also CONE, RAMP, SIMINTERF 
% 

%// $Revision: 1.3 $  $Date: 2001/09/28 14:24:32 $
%// Bert Kampes, 11-Dec-2000

%%% Handle input.
if     (nargin==2) ;
elseif (nargin==1) if (prod(size(nrows))==1) ncolumns=nrows;
		   else ncolumns=nrows(2); nrows=nrows(1); end;
elseif (nargin==0) nrows=512; ncolumns=nrows;
else   error('wrong number of input');
end;

%%% Pyramid generation (can smarter?)
x0      = round(ncolumns/2.);%          top of pyramid
y0      = round(nrows/2.);%             top of pyramid
[dx,dy] = meshgrid(1:ncolumns,1:nrows);
dx      = -abs(dx-x0);
dy      = -abs(dy-y0);
q       = find(dx>dy);
dx(q)   = dy(q);

%%% Normalize [-sqrt(x0):0] --> [0,1]
%minc  = min(dx(:));
minc   = -max(ncolumns-x0,nrows-y0);
dx     = (dx-minc) ./ -minc;

%%% EOF

