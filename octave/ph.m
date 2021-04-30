function c = ph(N)
%PH -- Shades of cyan and magenta color map.
%   PH(M) returns an M-by-3 matrix containing a hsv-like colormap.
%   It is simialr to the Stanford University 'mph' map,
%   and very suited for displaying interferometric phase.
%   (value==1, 3 revolving colors, 'transparant')
%   PH, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%     colormap(ph)
%
%   See also HSV, GRAY, COOL, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT, BRIGHTEN.

% Bert Kampes, 25/12/2000

% now create my own, brighten it... by beta=-0.2 (4x)
if (nargin<1) N = size(get(gcf,'colormap'),1); end
% for now only even, if odd, add one, subtract later...
if (rem(N,2)~=0) 
  warning('Sorry only even length cmap supported.');
  N=N+1;
end;

%%% hue
% difficult way.
%h = linspace(0,2.*pi,N);
%h = wrap(h)./(2.*pi) + 0.5;
% easy way
%h = linspace(.5,1.5,N+1).';
%h = h(1:N);
%assume even N, one double entry, periodic colors
h = linspace(.5,1,N/2);
h = [h, h-.5].';

%%% saturation
% easy way: interpolate continous vector. (normalized coordinates 0:1)
method = '*linear';
lo   = 0.35;% min. sateration
hi   = 0.6;%  max. sateration
%lo   = 0.5;% min. sateration
%hi   = 0.75;%  max. sateration
% how bout a cosine?
s_bp = [hi,lo,hi,lo,hi,lo,hi];% breakpoints
s_ax = (0:6)./6;
xi   = linspace(0,1,N);
s    = standing(interp1(s_ax,s_bp,xi,method));

%%%try a cosine for s, 3 peaks, bounded by hi,lo
%t=linspace(0,2*pi,N+1).';
%t=t(1:N);
%s = lo+(.5*(hi-lo)).*(cos(3*t) + 1);

%%% value
v = ones(N,1);

%%% return the cmap
c = hsv2rgb([h,s,v]);
%c = brighten(c,-0.2); %not very smart

%%% EOF

