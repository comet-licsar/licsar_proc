function [noise,gammageom,stdphase]=simnoise2(slope,Bperp)
% SIMNOISE  computes noise based on a matrix with slopes,
%    accounting for geometric decorrelation.
%    NOISE = SIMNOISE(SLOPE,BPERP) returns coherence for matrix
%    in radar coded system containing terrain slope and Bperp (scalar).
%    [NOISE,COH] = SIMNOISE(SLOPE,BPERP) returns coherence 
%    matrix as well.
%    [NOISE,COH,STDDEV] = SIMNOISE(SLOPE,BPERP) also returns matrix
%    with standard deviation.
%
%    This function is used in siminterf script.
%    Parameters used are valid for ERS1/2.
%

%// Erik Steenbergen, 22-Oct-2000
%// Bert Kampes, 07-Dec-2000
% Thanks goes to Ramon Hanssen for help.

%%% Set general parameters.
%% Fixed variables (do not change w/o reason).
Bw 	   = 15550000;%         [Hz] range bandwidth 
lambda     = 0.05666;%		[m] wavelength
theta      = deg2rad(20.);%	[rad] mean looking angle
R1         = 835000.;%		[m] mean range
c	   = 3e8;%	        [m/s] speed of light

% formulas courtesy of rh...
Bcritical = lambda*(Bw/c)*R1*tan(theta-slope);

% Due to trouble, adjust Bcritical matrix.
Bcritical(find(Bcritical<200))=200;
Bcritical(find(Bcritical>4000))=4000;
gammageom = abs(((Bcritical-Bperp)./Bcritical));
gammageom(isnan(gammageom)) = 0;
stdphase  = localphasestdev(gammageom);

%%% Compute gaussian noise, using correct sigma.
noise = randn(size(stdphase)) .* stdphase;

%%% EOF



%---------------------------------------------------
function [sdevphase] = localphasestdev(gam,N,figs);
% subfunction
% [sdevphase] = phasestdev(gam,[N,figs]);
%
%  INPUT
%    gam = gamma (abs coherence)
%    N   = number of effective looks
%    figs= 'yes'/'no' plot figures?
%  OUTPUT:
%    sdevphase [rad] standard deviation of the phase
%
% Computes the phase standard deviation for a given coherence and number
% of looks.  
%  If gamma is a vector or matrix, a lookup table is used.
%
% example, corresponding with 6-point cubic convolution as in
%          Hanssen and Bamler, 1999, TGARS
%     phasestdev(0.9988068,1,'yes');
%
% See also ...
%

% RH 16-Oct-2000 16:48 : Created from previous versions 
% RH 26-Oct-2000 15:45 : Added option for lookup table
%// BK 27-Oct-2000:      Changed lookup nearest to cubic, used as local.
%
if nargin == 0, help phasestdev;break;end
if nargin == 1, N=1; figs='y';end
if nargin == 2,      figs='y';end
if (prod(size(gam))<2) error('sorry only matrices'); end;

%%% matrix input
% BUILD LOOKUP TABLE
%gamma_table     = [0:0.1:0.9];
gamma_table     = [0:0.025:0.975];
sdevphase_table = zeros(size(gamma_table));
for ii = 1:length(gamma_table),
  points = 100;
  phi    = linspace(0,pi,points);
  dphi   = pi/(points-1);
  %    
  Y      = gamma_table(ii)*cos(phi);
  term1  = (1-gamma_table(ii)^2)^N/(2*pi);
  term2  = gamma(2*N-1)/( (gamma(N))^2*2^(2*N-2) );
  term3  = (2*N-1)*Y./(1-Y.^2).^(N+0.5).*acos(-Y)+(1-Y.^2).^(-N);
  term4  = 0;
  som    = 0;
  %    
  if N>1 
    term4 = 1/(2*N-2);
    for r = 0:N-2
      som = som+ gamma(N-0.5)*gamma(N-1-r)*(1+(2*r+1)*Y.^2) ./ ...
                 (gamma(N-0.5-r)*gamma(N-1)*(1-Y.^2).^(r+2));
    end
  end
  %   
  f2     = term1.*(term2.*term3+term4.*som);
  var2   = trapz(f2.*phi.^2);
  var2   = var2*dphi*2;
  sdevphase_table(ii) = sqrt(var2);
end

%%% Add the gamma==1, stdevphase==0 values
gamma_table     = [gamma_table, 1];
sdevphase_table = [sdevphase_table, 0];
% End building lookup table

%%% Start interpolation.
%method = 'nearest';
method = 'cubic';
for jj = 1:size(gam,1),
  sdevphase(jj,:) = interp1(gamma_table,sdevphase_table,gam(jj,:),method);
end


%%% EOF

