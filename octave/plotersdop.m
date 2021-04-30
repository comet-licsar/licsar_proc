function plotersdop(rawdatafile,len,nlines,headerbytes,ave);
% function plotersdop(rawdatafile,len,nlines,headerbytes,ave);
%
% Function  to estimate Doppler centroid using the average phase shift
%  from line to line
%
% Default values
%  ave = 15.5;
%  headerbytes = 412;
%  nlines = 400;
%  len    = 11644; % number of bytes
%

if nargin < 5,
  % ERS data are quantized in 5 bits [0:31]
  % Average = 15.5
  ave = 15.5;
  if nargin < 4,
    headerbytes = 412;
    if nargin < 3,
      nlines = 400;    % number of raw data lines
      if nargin < 2,
        len    = 11644; % number of bytes
        if nargin ==0,
          help plotersdop;
          break;
        end
      end
    end
  end
end
           
%rawdatafile = '25442.raw';
prf = 1679.902;
lambda = 0.0565;
vel    = 7554.267; %velocity m/s
theta  = 23 *pi/180; % look angle


datatype = 'int8';
%datatype = 'uint8';
fid = fopen(rawdatafile)       
fprintf(1,'Reading %i lines....\n',nlines);
%F=fread(fid,[nlines,11644],datatype);
F=fread(fid,[11644,nlines],datatype);
fclose(fid);

%Data start at byte 'headerbytes +1'
F = F(headerbytes+1:size(F,1),:);

% Assign complex values
dop = zeros(size(F,1)/2,1);

for n=1:2:nlines-1
signalreal = F(1:2:size(F,1),n) - ave;
signalimag = F(2:2:size(F,1),n) - ave;
signal1     = complex(signalreal,signalimag);
signalreal = F(1:2:size(F,1),n+1) - ave;
signalimag = F(2:2:size(F,1),n+1) - ave;
signal2     = complex(signalreal,signalimag);
dop = dop + signal2 .* conj(signal1);
end

mad = mean(angle(dop));
madprf = mean(angle(dop)/(2*pi));
 fprintf(1,'=> Average Doppler Centroid is %4.2f radians\n',mad);
 fprintf(1,'(2 pi radians correspond with 1 PRF)\n')
 fprintf(1,'(PRF ERS is %8.3f)\n',prf);
 fprintf(1,'=> That is %4.2f of the prf\n', madprf);
 fprintf(1,'=> Estimated Doppler Centroid Frequency %4.2f Hz\n',madprf*prf);

% Squint angle
 squint = asin( ((madprf*prf) * lambda)/(2*vel*sin(theta)));
 fprintf(1,'=> The squint angle is %4.2f degrees\n',squint*180/pi);

figure(1);plot(angle(dop)/(2*pi));
xlabel('Range point (bin) number');
ylabel('Doppler Centroid frequency expressed as fraction of the PRF')
