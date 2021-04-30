function [data,header,cmap] = freadras(filename);
% FREADRAS  --  Read a SUNraster file.
%
% [DATA, HEADER, CMAP] = FREADRAS(filename);
%
% FREADRAS reads 8 and 32 bit depth rasterf files.
%

%// Bert Kampes, 06-Oct-2000
more off

%%% Check input.
if ( nargin == 0 ) help readras; error('not enough input'); end;
if (~ischar(filename))
  error('freadras: argument must be a string');
end;
fid = fopen(filename);
if (fid<0) error(ferror(fid)); end;

%%% The header
disp(['Reading sunraster header 32B for file: ',filename]);
header = fread(fid,[8 1],'int32');
disp(['ras_magic:     ', int2str(header(1))]);
disp(['ras_width:     ', int2str(header(2))]);
disp(['ras_height:    ', int2str(header(3))]);
disp(['ras_depth:     ', int2str(header(4))]);
disp(['ras_length:    ', int2str(header(5))]);
disp(['ras_type:      ', int2str(header(6))]);
disp(['ras_maptype:   ', int2str(header(7))]);
disp(['ras_maplength: ', int2str(header(8))]);

%%% The cmap.
if (header(8)~=0)
  disp(['Reading cmap from raster file']);
  r    = fread(fid,[header(8)/3 1],'uint8');
  g    = fread(fid,[header(8)/3 1],'uint8');
  b    = fread(fid,[header(8)/3 1],'uint8');
  cmap = [r g b];
  cmap = cmap ./ 256;
else
  disp(['No cmap in raster file data']);
end

%%% The data.
disp(['Reading sunraster data']);
switch header(4)
  case 8
    data = fread(fid,[header(2) header(3)],'uint8');
  case 32
    data = fread(fid,[header(2) header(3)],'float32');
  otherwise
    error(['Unknown format, data depth = ',num2str(header(4))]);
end;


fclose(fid);
% major col order in matlab, so transpose
data=data.';


%%%EOF
more on

