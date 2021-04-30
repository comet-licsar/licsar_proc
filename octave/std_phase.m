function std_phase = std_phase (coherence,L);
% STD_PHASE -- Phase standard deviation for single look distributed scatterers.
%
%   STD_PHASE(COHERENCE) returns the standard deviation for the phase
%   based on the coherence. This should be valid for distributed 
%   scatterers. COHERENCE is a (complex) number with absolute value
%   between 0 and 1. If COHERENCE is a vector, a vector is returned.
%   Note that this is not the standard deviation of a normal distribution.
%
%   STD_PHASE(COHERENCE,L) returns std for multilookfactor L. (NOT YET)
%
%   Output is in degrees.
%
%   General equation for std is (expection E{.}):
%     var = int_{-pi}^{+pi} [phi-E{phi}].^2 .* pdf(phi) dphi
%   For single look this can be written in closed form as (coherence G):
%     var = 1/3*pi^2 - pi*arcsin(G) + arcsin^2(G) - 0.5*Li(G^2)
%   where Li(G)=sum_1^infty G^2k./k^2 (Euler's dilogarithm). 
%
%   Example:
%     std_phase(0.8);% return std of phase for coherence of 0.8    
%   Or a nice plot:
%     coh=0:0.01:1; plot(coh,std_phase(coh));xlabel('coherence');ylabel('std');
%     title('phase standard deviation vs. coherence'); grid on;
%
%   See also: PDF_PHASE, 
%

%   Or for multiple levels of multilooking:
%     coh=0:0.01:1;
%     plot(coh,std_phase(coh,1),'k',coh,std_phase(coh,10),'b',...
%          coh,std_phase(coh,20),'r');
%     xlabel('coherence');ylabel('std');
%     title('phase standard deviation vs. coherence'); grid on;
%

%     (Hanssen 2001, eq. (4.2.28), p. 94)
%

% Rens Swart * 17 april 2001
%// BK 18-Apr-2001
%// $Revision: 1.2 $  $Date: 2001/09/28 14:24:33 $
%
% Thinking of renaming this function to coh2std.


%%% Check input.
exitwithhelp=0;
switch nargin
  case 1
    L=1;
  case 2
    disp('sorry not yet');
  otherwise
    coherence=0.5; L=1;% (dummies)
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
if (exitwithhelp==1) helphelp; break; end;


%%% Compute std.
switch (L)
  case 1
    var_phase = (pi.^2)./3 - pi.*asin(coherence) + ...
                asin(coherence).^2 - Li2(coherence)./2;
  otherwise
    error('only single look for now, sorry')
end;

%%% what to return?
std_phase = rad2deg(sqrt(var_phase));

%%% EOF.



%---------------------------------------------------------------
% SUBFUNCTION LI2
%---------------------------------------------------------------
function Li2 = Li2 (coherence)
% Li2 = Li2 (coherence, kmax)
%   Calculate Euler's dilogarithm for use in calculation of
%   phase variance.
%   Note that actually the argument of Euler's dilogarithm should be
%   the squared absolute coherence. This function however has the
%   absolute coherence as parameter.
%   (Hanssen 2001, eq. (4.2.29), p. 95)
%   Li(G)=sum_1^infty G^2k./k^2 (Euler's dilogarithm). 
%
%   Implementation with loop cut of at 100, which is ok, 
%   even for coherence=0.99 (BK).
%   If input is a vector, slow implemented? vector returned.
%   coherence should be lying...
%

% Rens Swart * 17 april 2001
%if (nargin~=1) error('only 1 argument in subfunction Li2.'); end;
%if (coherence<0 | coherence>1) error('coherence not in [0,1]'); end;
%%% check input is vector/scalar

coherence = lying(coherence);
KMAX = 100;
k    = standing(1:KMAX);
G    = (ones(KMAX,1)*coherence).^(2.*k*ones(1,length(coherence)));
K    = k.^2*ones(1,length(coherence));
Li2  = sum(G./K,1);

%%% TEST
%for ii=1:length(coherence)
%  g    = coherence(ii);
%  Li2a = sum((g.^(2.*k)) ./ (k.^2))
%end
%%% EOF.


