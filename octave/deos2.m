function c = deos2(N);
% DEOS2 --  colormap for wrapped interferometric phase.
%   DEOS2(M) returns an M-by-3 matrix containing a blue-green-red-blue colormap.
%
%   DEOS2, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%     colormap(deos2)
% 
%   See also HSV, GRAY, COOL, HOT, BONE, COPPER, PINK, FLAG, PH, DEOS,
%   COLORMAP, RGBPLOT, BRIGHTEN.
 
%// $Revision: 1.1 $  $Date: 2001/03/16 13:47:37 $
%// Bert Kampes, 27-Dec-2000

%%% Handle input.
if (nargin<1) N=size(get(gcf,'colormap'),1); end

%%% rgb breakpoints.
bp = [0,   1/3    2/3   1];
r  = [0,   1,     0,    0];
g  = [0,   0,     1,    0];
b  = [1,   0,     0,    1];
 
%%% Interpolate colormap to new size.
np = linspace(0,1,N);
r2 = interp1(bp,r,np,'linear');
g2 = interp1(bp,g,np,'linear');
b2 = interp1(bp,b,np,'linear');
 
%%% Return colormap.
c = [r2.',g2.',b2.'];
c = brighten(c,0.6);
%rgbplot(c);

%%% EOF

