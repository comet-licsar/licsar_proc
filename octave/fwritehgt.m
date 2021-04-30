function count = fwritehgt(ampli,phase, outfile);
% FWRITEHGT writes amplitude and phase to binary hgt file.
%
%   COUNT = FWRITEHGT(ampli, phase, filename) writes amplitude, phase matrix
%   to specified file in hgt format (band interleaved, major row order).
%   returns number of elements succesfully written. (?)
%
%   See also FREADBK, FWRITEBK, FREADHGT, FREAD, FOPEN, FWRITE
%

% $Revision: 1.4 $  $Date: 2000/07/19 14:22:08 $
% Bert Kampes, 7/3/00

%%% Handle input.
if (nargin <  2) error('FWRITEHGT: no data specified.'); end;
if (nargin <  3)
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
end;
if (size(ampli) ~= size(phase))
  error('WRITEHGTFILE: size ampli not equal to size phase.'); end;

%%% Write data to file in major row order.
% Could use here: count = fwritebk([ampli,phase], outfile);
% but mem?
%
fid = fopen(outfile,'w');
if (fid<0)%				try one more time.
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
  fid = fopen(outfile,'w');
  if (fid<0) error('FWRITEHGT: outfile could not be opened.'); end;
end;

% Written sequentially (?)
count=0;
[lines, width] = size(ampli);
for ii=1:lines
  cnt=fwrite(fid,ampli(ii,:),'float32');
  cnt=fwrite(fid,phase(ii,:),'float32');
  count=count+cnt;
end
fclose(fid);

%%% EOF








