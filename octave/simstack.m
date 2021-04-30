function [DATACUBE,Q,V,trend] = simstack(B,T,Q,V,trend,noise);
%  SIMSTACK  -- Simulate a datacube with phase (use approximations).
%
%
%  Input:
%    DATACUBE = SIMSTACK(B,T) simulates phase for K interferograms
%      with baselines B,T size (azi,range)=(200x250).
%      DATACUBE(ii,:,:) is phase for IFG ii (azimuth,range).
%
%    SIMSTACK(b,t,Q) uses Q as matrix for height; Q is DEM in meters in
%      radarcoordinates (azimuth,range).
%       Actually better to see it as residual topography error, so scaled to [-Q,Q].
%       (i.e., (4pi/lambda)*(B/r*sin(theta))*Q is used to make topogr. phase.)  
%       An incidence angle of 23 degrees is used, range=850000.
%      if Q is a scalar then fractal is used with heights [-Q,Q].
%      if Q==0 then ignored;
%      default Q=10.
%
%    SIMSTACK(b,t,q,V) uses V as unit matrix for deformation; V is line
%      of sight linear velocity rate matrix in mm (for one year if T is in years).
%      (i.e., 4pi/lambda*T.*V is used to generate velocities phase.)
%      If V is a scalar then fractal is used with velocities [0,V].
%      If V==0 then ignored.
%      default V = -20. (subsidence)
%
%    SIMSTACK(b,t,q,v,TREND) uses TREND to add trends to ifgs.
%      If TREND is a matrix (Kx3) then this is used to add a plane with 
%       bias TREND(ii,1), and TREND(ii,2) fringes in azimuth and TREND(ii,3)
%       fringes in range.
%      Such a TREND matrix can also be returned (if scalar input, see below).
%      If TREND is a scalar, then a trend is added to each IFG with the
%       specified std.dev. in number of fringes, and a random bias.
%      If TREND==0 then ignored.
%      Default no TREND.
%
%    SIMSTACK(b,t,q,v,trend,NOISE) uses NOISE as stddev of gaussian noise.
%      Noise in degrees. added is randn(size(Q));
%      If NOISE==0 then ignored.
%      Default no NOISE.
%
%
%  Output:
%    [DATA, Q, V, trend] = ... returns matrices Q,V,trend as well.
%
%
%  Examples:
%    To create a stack of 6 interferograms with maximum DEM error of 5 meters,
%    subsidences of maximum 10 mm/year, no trends and noise:
%      M  = 6;  % number of acquisitions
%      dB = 500;% sigma of baselines distr. (expect 95% within 2sigma)
%      dT = 5;  % years of acquisitions
%      T  = sort(rand(M,1)*dT);
%      B  = randn(M,1)*dB; q=find(abs(B)>1400); B(q)=randn(size(q))*(dB/2);
%      %plot(B,T,'r+'); title('IFG distribution');
%      D  = simstack(B,T,5,-10);
%      imagesc(squeeze(D(1,:,:))); colorbar; title('IFG1');
%      
%    With the same B,T distribution, but with noise, use
%      noise = 20.;% degrees
%      D  = simstack(B,T,5,10,0,noise);
%
%  See also:
%    DEOS InSAR and fractal toolboxes
%

%%% TODO: atmo fractal based, check validity of formulas, test
%%%  (may be useful to simulate [-V:+V], no trends)
%%%  since in ps_solveAPS an APS is estimated that covers this?
%%%  This already can be done by creating an input V...

%// BK 15-May-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:46 $

%%% Check input
n_azi   = 200;%		[pixels] default, 800m single look.
n_range = 250;%		[pixels] default, 5000m single look.
if (nargin<2 | nargin>6) error('wrong #of input'); end;
if (nargin<6)
  noise=0;
end;
if (nargin<5)
  trend=0;
end;
if (nargin<4)
  V=-20;%				mm/year in LOS.
end;
if (nargin<3)
  Q=10;%				default fractal with 10 m. DEM errors
end;
%%%
if (~exist('fracsurf')) error('fracsurf not known, download fractal toolbox.'); end;
if (~exist('ramp')) error('ramp not known, download insar toolbox.'); end;
%%%
K         = length(B);%				number of IFGs
if (length(B)~=length(T)) error('length B,T not equal'); end;

%%% Check scalar input to generate fractal surfaces.
%%% Check matrix input.
if (prod(size(Q))~=1)
  [n_azi,n_range]=size(Q);
end;
if (prod(size(V))~=1)
  [n_azi,n_range]=size(V);
end;
if (prod(size(Q))~=1 & prod(size(V))~=1)
  if (any(size(Q)~=size(V)))
    n_azi   = min(size(Q,1),size(V,1));
    n_range = min(size(Q,2),size(V,2));
    Q = Q(1:n_azi,1:n_range);
    V = V(1:n_azi,1:n_range);
    warning(['size input matrix Q and V not same, using smallest: ',...
	    num2str([n_azi,n_range])]);

  end
end
%%% Generate topo/defo with fractals.
if (prod(size(Q))==1)
  Q       = (2.*Q).*(topo(n_azi,n_range) - 0.5);
end
if (prod(size(V))==1)
  V       = V.*defo(n_azi,n_range);
end
if (prod(size(trend))==1)
  if (trend~=0)
    trend = trend.*randn(K,2);%			nfringes in azi/range per IFG
    trend = [(rand(K,1)-.5).*2.*pi trend];%	bias [-pi,pi]
  else 
    trend = zeros(K,3);
  end
end


%%% Local variables.
[n_azi, n_range] = size(Q);
lambda    = 0.0566;
pi4       = 4*pi;
alpha     = deg2rad (23);%			local incidence angle
rsintheta = 850000*sin(alpha);
DATACUBE  = zeros(K,n_azi,n_range);



%%% Simulate phase per IFG.
for ii=1:K
  DATACUBE(ii,:,:) = wrap(...
		     B(ii).*(pi4./(lambda.*rsintheta)).*Q + ...	%topographic phase
		     T(ii).*(pi4./lambda).*1e-3.*V        + ...	%deformation phase
		     trend(ii,1)                          + ... %bias
		     trend(ii,2).*2.*pi.*ramp(n_range,n_azi).'   + ... %fringes in azi
		     trend(ii,3).*2.*pi.*ramp(n_azi,n_range)     + ... %fringes in range
		     deg2rad(noise).*randn(size(Q)));		%gaussian noise
end
if (K==1) DATACUBE=squeeze(DATACUBE); end;
%%% EOF


%-----------------------------------------------------
% subfunction generate topo matrix (scaled to [0,1])
%-----------------------------------------------------
function T = topo(n_azi, n_range, beta)
if (nargin>3)  error('wrong input'); end;
if (nargin<3)  beta    = 1;   end;%	default little covariance in DEM errors
if (nargin<2)  n_range = 100; end;
if (nargin==0) n_azi   = 100; end;
%
T = fracsurf(max(n_azi,n_range), beta ,'n',100);
T = T(1:n_azi,1:n_range);
T = T-min(T(:));%			[0,x]
T = T./max(T(:));%			[0,1]
%%% EOF

%-----------------------------------------------------
% subfunction generate defo matrix per unit of time (scaled to [0,1])
%-----------------------------------------------------
function D = defo(n_azi, n_range, beta)
if (nargin>3)  error('wrong input'); end;
if (nargin<3)  beta    = 4;   end;%	default smooth defo pattern.
if (nargin<2)  n_range = 100; end;
if (nargin==0) n_azi   = 100; end;
%
D = fracsurf(max(n_azi,n_range), beta ,'n',100);
D = D(1:n_azi,1:n_range);
%D = D/ -min(D(:));%			normalize [-1,0]
%D(find(D>0))=0;%			no uplift
D = D-min(D(:));%				[0,x]
D = D./max(D(:));%			[0,1]
%%% EOF


