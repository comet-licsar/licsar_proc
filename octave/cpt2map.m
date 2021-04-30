function cmap = cpt2map(cptfile,numcolors);
% CPT2MAP -- convert a cpt table to a matlab colormap.
%  
%   A cpt table is used in GMT. It contains entries which couple 
%   a data range to RGB values [0:255]. Matlab colormaps contain 3 
%   columns with RGB values in the interval [0:1].
%   The matlab colormap is uniformly spaced, linearly interpolated.
%   .cpt is the default extension, which can be omited.  
%  
%   CMAP = CPT2MAP(CPTFILE) returns map of length 64.
%   CMAP = CPT2MAP(CPTFILE,SIZE) return map of length SIZE.
%  
%   Cpt format is defined as ascii files with lines:
%     z0    Rmin  Gmin  Bmin   z1     Rmax  Gmax  Bmax [A]
%   For example:
%     ...
%     0.17   143   143   192   0.54   159   159   124
%     0.54   159   159   124   0.91   175   175    56
%     ...
%
% See also BRIGHTEN, COLORMAP, http://imina.soest.hawaii.edu/gmt/
% INTERP1
%

%// $Revision: 1.2 $  $Date: 2001/03/16 13:46:36 $ 
%// Bert Kampes, 22-Dec-2000


%%% Handle input options.
if (nargin==1)
  numcolors=64;% default
elseif (nargin==2)
  ; %do nothing
else
  warning('no arguments');
  help cpt2mat;
  break;
end
if (~ischar(cptfile)) error('input argument 1 has to be a string'); end;

%%% Obtain cptfile, try filename .cpt
if (~exist(cptfile,'file')) cptfile=[cptfile,'.cpt']; end;
if (~exist(cptfile,'file')) 
  warning(['cpt file ', cptfile, ' not found. Please select file.']);
  [infile, inpath] = uigetfile('*', 'Select cpt file', 0,0);
  cptfile = [infile, inpath];
end;
cmap = load(cptfile);

%%% Translate cpt to cmap, remove extension (better mean of 1,5 and 2:4,6:8 ?)
xaxis = cmap(:,1);%   z0 data entry
cmap  = cmap(:,2:4);% RGB entries

%%% Rescale [0:255] --> [0:1]
cmap = cmap ./ 255;
%if (numcolors==0) numcolors=size(cmap,1); end

% Interpolate to uniform grid and optionally to new size.
method = 'linear';
xi     = linspace(xaxis(1),xaxis(length(xaxis)),numcolors);
cmap   = interp1(xaxis,cmap,xi);

%%% EOF

