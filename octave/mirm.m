function out = mirm(varargin);
% MIRM  --  Multi-Image Reflectivity Map.
%
%   M = MIRM(datacube) returns a matrix containing the mean amplitude
%   of the interferograms in the datacube.  Each plane of the cube 
%   contains an ifg.  CUBE(ii,ll,pp) where ii is the ifg number, ll the 
%   azimuth, and pp the range coordinate.
%   M = MIRM(datacube,FACTOR) oversamples in range by factor.
%
%   M = MIRM(matrix1, matrix2, ...) returns a matrix containing the mean
%   magnitude.
%
%   M = MIRM(matrix1, matrix2, ..., FACTOR) optionally uses the scalar
%   FACTOR to oversample the data in the Y direction.
%
%   If data is complex, magnitude is taken before anything else. (correct?)
%   If oversampling is requested, this is done as the last step. (correct?)
%   FFT method is used for interpolation (oversampling).
%
%   Note that calibration of the SLC images is required in general.
%
%   See also INSAR toolbox.
%
%
% WARNING: it seems to be better to oversample the COMPLEX ifgs first,
% before computing the average.  if you do interpolation as last step,
% then a negative amplitude may result!
% but this would also happen if you offer magnitude images to this function
% so what is the difference?
%

%// BK 07-Aug-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:45 $

numargs = length(varargin);
if (numargs==0 | isscalar(varargin{1})) helphelp; break; end;


%%% Check for FACTOR as last argument
FACTOR = 1;
if (isscalar(varargin{numargs}))
  FACTOR = varargin{numargs};
  numargs = numargs - 1;
end;


%%% Check if input was datacube.
dc = 0;% no datacube as input
%if (nargin==1 & size(varargin{1},3)~=1) dc = 1; end;
%if (nargin==2 & size(varargin{1},3)~=1) dc = 1; end;
if (size(varargin{1},3)~=1) dc = 1; end;% assume 3D datacube now

%%% DEBUG
%numargs
%dc
%FACTOR
%keyboard

FIRSTOVERSAMPLE=1;
if (FIRSTOVERSAMPLE==1)
  disp('first oversampling complex images in range');
  warning('only for datacube input, testing');
  rangeindex = 3;% for datacube, 3rd dim. is range.
  %out = interpft(out,FACTOR*size(out,rangeindex),rangeindex);
  %varargin{1} = interpft(varargin{1},FACTOR*size(varargin{1},rangeindex),rangeindex);
  %%% seems to be a problem with interpft on arrays? losing a dimension? so:
  out = varargin{1};
  out2 = zeros(size(varargin{1},1),size(varargin{1},2),FACTOR*size(varargin{1},3));
  %keyboard
  for ii=1:numargs
    q = squeeze(out(ii,:,:));
    rangeindex = 2;% squeezed to 2 dimensions.
    out2(ii,:,:) = interpft(q,FACTOR*size(q,rangeindex),rangeindex);
  end
  varargin{1} = out2; clear out out2;
end
%%% Compute mean amplitude.
if (~isreal(varargin{1}))
  varargin{1} = abs(varargin{1});
end;
if (dc==1)
  ifgindex = 1;% 1 is dimension over which ifgs are.
  out = squeeze(mean(varargin{1},ifgindex));

%%% Lot of matrices as input.
else
  out = varargin{1};
  for ii = 2:numargs;
    if (isreal(varargin{ii}))
      out = out + varargin{ii};
    else
      out = out + abs(varargin{ii});
    end
  end
  out = out ./ numargs;
end


if (FIRSTOVERSAMPLE~=1)
%%% Oversample result. (applied to each range; dimension 2)
%%% Use FFT method.
rangeindex = 2;% for mean map, which is squeezed, 2nd dim is range.
out = interpft(out,FACTOR*size(out,rangeindex),rangeindex);

%%% Make sure after interpolation no values are below zero...?
q = find(out<0);
if (~isempty(q))
  warning('after interpolation some values<0 (i set them to 0)');
  out(q)=0;
end
end

%%% PLOT temporarly...
figure
scale=150;
q=out.^0.3;
q=scale.*q./mean(q(:));
q(find(q>255)) = 255;
q(find(q<16))  = 16;
imagesc(q,[16 255]);
colormap(gray);
colorbar

%%% EOF

