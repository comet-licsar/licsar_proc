function std_phasePS = std_phasePS (coherence);
% STD_PHASEPS -- Phase standard deviation for point scatterers.
%
%   STD_PHASEPS(COHERENCE) returns the standard deviation for the phase
%   based on the coherence. This should be valid for distributed 
%   scatterers. COHERENCE is a (complex) number with absolute value
%   between 0 and 1. If COHERENCE is a vector, a vector is returned.
%   Note that this is not the standard deviation of a normal distribution.
%
%   STD_PHASEPS(COHERENCE,L) returns std for multilookfactor L.
%   (NOT possible for ponit scatterers...)
%
%   Output is in degrees.
%
%   For point scatterers the variance of the phase is (Rodriguez&Martin 1992)
%     var = (1-g^2)./(2g^2);
%
%   Example:
%     std_phasePS(0.8);% return std of phase for coherence of 0.8    
%
%   See also: PDF_PHASE, STD_PHASE 
%

%     (Hanssen 2001, eq. (4.2.28), p. 94)
%

%// BK 18-Apr-2001
%// $Revision: 1.2 $  $Date: 2001/09/28 14:24:33 $
%


%%% Check input.
exitwithhelp=0;
switch nargin
  case 1
    L=1;
  case 2
    warning('sorry no multilooking for point scatterers');
  otherwise
    coherence=0.9; L=1;% (dummies)
    exitwithhelp=1;
end
if (~isreal(coherence)) coherence=abs(coherence); end;
if (min(size(coherence))~=1)
  warning('coherence scalar or vector.');
  exitwithhelp=1;
end;
if (min(coherence)<0 | max(coherence)>1) 
  warning('coherence not in [0,1]');
  exitwithhelp=1;
end;
if (coherence<0.8) warning('a point scatterer with low coherence?'); end;
if (exitwithhelp==1) helphelp; break; end;


%%% Compute std.
%var_phase = (1-coherence.^2)./(2.*L.*coherence.^2);
var_phase = (1-coherence.^2)./(2.*coherence.^2);

%%% what to return?
std_phasePS = rad2deg(sqrt(var_phase));

%%% EOF.

