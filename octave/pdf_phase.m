function pdf = pdf_phase(phase,coherence,L,phase_0)
%PDF_PHASE  --  Probability density function for single look phase data.
%
%   Probability density function (pdf) for the phase is expressed    
%   as function of phase, coherence, and expected center.
%   PDF_PHASE(PHASE) returns the pdf for COHERENCE=0.9, PHASE_0=0.
%     PHASE input is in radians.
%   PDF_PHASE(PHASE,COHERENCE) uses the specified coherence.
%   PDF_PHASE(PHASE,COHERENCE,L) uses the specified multilookfactor.
%   PDF_PHASE(PHASE,COHERENCE,L,PHASE_0) uses the specified center.
% 
%   PHASE should be scalar or lying vector (x-axis).
%   COHERENCE can be complex, 0<=abs(coh)<1. 
%   COHERENCE and PHASE_0 either scalars or standing vectors,
%   but most be same size.
%   Valid for distributed scattering.
%
%   Output is a matrix with rows for coherence, columns over phase.
%
%   Based on formulas in thesis R. Hanssen.
%   (References to Just and Bamler, 1994; Tough et al., 1995):
%   For L==1:
%     pdf(p,g,p0) = (1-abs(g)^2)./2pi *
%                 1./(1-abs(g).^2*(cos(p-p0)).^2) *
%                 (1+((abs(g)cos(p-p0)arccos(-abs(g)*cos(p-p0))) ./ 
%                     (sqrt(1-abs(g).^2*(cos(p-p0)).^2)))))
%   For L~=1:
%     pdf(p,g,L,p0) = 1/2pi*(1-g^2)^L {...etc., see source}
%
%   Example to view a shifted multi-modal pdf for single look data:
%     phase=-pi:0.01:6*pi;
%     coher=0.7;
%     pha_0=deg2rad(60);
%     pdf=pdf_phase(phase,coher,1,pha_0);
%     plot(rad2deg(phase),pdf);
%     title(['pdf for coherence=',num2str(coher)]);
%     xlabel('phase [deg]');
%     ylabel('pdf(phase)');
%
%   Compare multilooking effect for same coherence:
%     phase = linspace(-pi,pi);
%     coher = 0.7;
%     pdf01=pdf_phase(phase,coher,1);
%     pdf20=pdf_phase(phase,coher,20);
%     plot(phase,pdf01,'b',phase,pdf20,'r');
%     title('pdf for coherence=0.7; single look (blue); multilook L=20 (red)');
%
%   Example to plot pdf for number of coherence levels:
%     phase=linspace(-pi,pi);
%     coher=standing(0.1:0.2:0.9);
%     pdf=pdf_phase(phase,coher);
%     plot(phase,pdf);
%     legend(num2str(coher)); title('single look pdfs for coherence levels'); 
%     grid on;
%

%// BK 19-Feb-2001
%// $Revision: 1.3 $  $Date: 2001/09/28 14:24:32 $
%
% Thinking of renaming this function to coh2pdf.

%%% Handle input.
exitwithhelp=0;
if (nargin>4) warning('too much input');exitwithhelp=1; end;
if (nargin<1) warning('not enough input');phase=1;exitwithhelp=1; end;
if (nargin<2) coherence=0.9; end;
if (nargin<3) L=1; end;
if (nargin<4) phase_0=zeros(size(coherence)); end;
if (max(abs(coherence))>=1)
  warning('pdf_phase: abs(coherence)>=1');
  exitwithhelp=1;
end;
%
% Check scalar/vector.
%if (prod(size(phase))~=max(size(phase)))
if (size(phase,1)~=1)
  warning('pdf_phase: input1 phase should be scalar or lying vector.');
  exitwithhelp=1;
end;
if (size(coherence,2)~=1)
  warning('pdf_phase: input2 coherence should be scalar or standing vector.');
  exitwithhelp=1;
end;
if (~isequal(size(coherence),size(phase_0)))% same dimensions?
  warning('pdf_phase: input3 phase_0 should be same size as input2 coherence.');
  exitwithhelp=1;
end;
if (L<1 | rem(L,1)~=0)
  warning('3rd argument integer multilookfactor');
  exitwithhelp=1;
end;
if (exitwithhelp==1) helphelp; break; end;


%%% Compute pdf.
%g     = abs(coherence);% could be complex.
%gcosp = g*ones(1,length(phase)).* ...
%        cos((ones(length(g),1)*phase)-(phase_0*ones(1,length(phase))));
%g2    = g*ones(1,length(phase));
g2    = abs(coherence)*ones(1,length(phase));
gcosp = g2.*cos((ones(length(coherence),1)*phase)-(phase_0*ones(1,length(phase))));

switch L
  case 1
    %%% input with vectors of coherence/phase_0, does not seem efficient...
    %pdf   = ((1-g.^2)*ones(1,length(phase))) ./ ...
    pdf   = (1-g2.^2) ./ ...
            ((2.*pi).*(1-(gcosp).^2)) .* ...
            (1+((gcosp.*acos(-gcosp)) ./ (sqrt(1-gcosp.^2))));
  otherwise
    part1 = ((1-g2.^2).^L)./(2.*pi);
    part2 = gf(2.*L-1)./(gf(L).^2*2.^(2.*L-2));
    part3 = (((2.*L-1)*gcosp)./((1-gcosp.^2).^(L+0.5))) .* ...
	      (pi./2+asin(gcosp)) + ...
              (1-gcosp.^2).^-L;
    part4 = comppart4(L,gcosp);
    pdf   = part1.*(part2*part3+part4);
end

%%% EOF.


%----------------------
% SUBFUNCTION
%----------------------
function p = comppart4(L,gcosp);
% part of total, see tough et al.
r = 0:L-2;
s = zeros(size(gcosp));
for r=0:L-2
  s = s + ( (gf(L-0.5)./gf(L-0.5-r)) .* ...
	  (gf(L-1-r)./gf(L-1)) .* ...
	  ((1+(2.*r+1)*gcosp.^2)./((1-gcosp.^2).^(r+2))) );
end
p = (1./(2.*(L-1))) .* s;
%%% EOF.



%----------------------
% SUBFUNCTION
%----------------------
function G = gf(L)
%gamma(x) = integral from 0 to inf of t^(x-1) exp(-t) dt.
% gamma function for integers: gf(L)=(L-1)!
%G = prod(1:L-1);
G=gamma(L);% exists already in matlab...
%%% EOF.



