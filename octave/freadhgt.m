function [phase,ampli] = freadhgt(infile, numlines);
% [PHASE, AMPLI] = FREADHGT(HGTFILE, NUMLINES);
%   Read specified HGTFILE into matrices for PHASE and AMPLItude.
%

% this will be a format option in freadbk, as will 'mph' in future
% it simply reads in band interleaved?
% [phase, ampli] = ...   to prevent ampli return if not requested.
%

% $Revision: 1.3 $  $Date: 2000/05/03 11:13:50 $
% Bert Kampes, 08-Mar-2000



%%% Handle input
if (nargin == 0)
  [infile, inpath] = uigetfile('*', 'Select Data File', 0,0);
  numlines = input('Enter number of lines: ');
  infile = [infile, inpath];
elseif (nargin ==1)
  numlines = input('Enter number of lines: ');
end


%%% Open file, better line by line for memory?
phase = freadbk(infile,numlines,'float32');
[L P] = size(phase);
if (rem(P,2) ~= 0) error('size hgt file not ok'); end;
if (nargout==2) ampli = phase(:,1:P/2); end;
phase = phase(:,P/2+1:P);

%%% EOF
