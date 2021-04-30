function slc = simulateslc(lines, pixels);
% SIMULATESLC  generate a SLC image
%   SIMULATESLC(N) returns a simulated Single Look Complex 
%   radar image of size NxN.
%   SIMULATESLC(N,M) returns size NxM.
%
%   phase is simulated as: 2*pi*rand(N,M)
%   ampli is simulated as: sqrt(-log(rand(N,M)))
%
%   not band limited or shifted
%   See also RAND, FRACSURF,
%

% $Revision: 1.3 $  $Date: 2001/03/16 13:47:05 $
%// Bert Kampes, 16-Jun-2000

%%% Handle input
if ( argv == 1 ) pixels = lines; end;

% simulate data rule of paper
% Ramon Hanssen and Richard Bamler.
% Evaluation of interpolation kernels for SAR interferometry.
% IEEE Trans. on Geoscience and Remote Sensing,
% 37(1):318-321, January 1999.
%phase = 2*pi*rand(lines,pixels);
%ampli = sqrt(-log(rand(lines,pixels)));
%slc   = ampli .* exp(i*phase);
slc   = sqrt(-log(rand(lines,pixels))) .* exp(i*2*pi*rand(lines,pixels));
%DATA  = fftshift(fft(data));


%%% EOF
