function c = ramp(nrows,ncolumns)
% RAMP --  Generate unit ramp (left to right).
%   c = RAMP;        ramp [0,1] of dimensions (512,512);
%   c = RAMP([L P]); ramp [0,1] of dimensions (L,P);
%   c = RAMP(L);     ramp [0,1] of dimensions (L,L);
%   c = RAMP(L,P);   ramp [0,1] of dimensions (L,P);
% 
%   To create a wrapped phase ramp with 4 wraps, use
%     c = wrap(4*2*pi*ramp(512,128)); imagesc(c); colorbar;
%
%   To create a phase trend for an interferogram with 4 fringes in X
%   and 2 fringes in Y, one could use:
%     Y  = 512; X  = 128;
%     nY = 2;   nX = 4;
%     cX = nX*2*pi*ramp(Y,X);
%     cY = nY*2*pi*ramp(X,Y);
%     c  = wrap(cX + cY.'); imagesc(c); colorbar;
%
%   See also PYRAMID, CONE, SIMINTERF 
% 

%// $Revision: 1.2 $  $Date: 2001/09/28 14:24:32 $
%// Bert Kampes, 11-Dec-2000

%%% Handle input.
%ramp (nrows, ncolumns)
%        1      2
if     (nargin==2) ;
elseif (nargin==1) if (prod(size(nrows))==1) ncolumns=nrows;
		   else ncolumns=nrows(2); nrows=nrows(1); end;
elseif (nargin==0) nrows=512; ncolumns=nrows;
else   error('wrong number of input');
end;

%%% Ramp generation (to easy for function?)
maxheight = 1;
c         = ones(nrows,1) * lying(linspace(0,maxheight,ncolumns));

%%% EOF

