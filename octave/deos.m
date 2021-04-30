function c = deos(N);
% DEOS --  colormap for wrapped interferometric phase.
%   DEOS(M) returns an M-by-3 matrix containing a hsv-like colormap.
%
%   DEOS, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%     colormap(deos)
% 
%   See also HSV, GRAY, COOL, HOT, BONE, COPPER, PINK, FLAG, PH,
%   COLORMAP, RGBPLOT, BRIGHTEN.
 
%// $Revision: 1.2 $  $Date: 2001/03/16 13:46:39 $
%// Bert Kampes, 27-Dec-2000

%%% Handle input.
if (nargin<1) N=size(get(gcf,'colormap'),1); end

%%% rgb breakpoints.
bp = [0,   0.25,  0.5,  0.75,  1];
r  = [1,   0,     0.5,  1,     1];
g  = [1,   1,     0,    0.5,   1];
b  = [0,   1,     1,    0.5,   0];
 
%%% Interpolate colormap to new size.
np = linspace(0,1,N);
r2 = interp1(bp,r,np,'linear');
g2 = interp1(bp,g,np,'linear');
b2 = interp1(bp,b,np,'linear');
 
%%% Return colormap.
c = [r2.',g2.',b2.'];
hs = rgb2hsv(c);
hs(:,3)=1; %value==1
c = hsv2rgb(hs);

%%% EOF

