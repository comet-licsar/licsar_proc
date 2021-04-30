function ha = heightamb(bperp)
% HEIGHTAMB -- Approximation of height ambiguity for perpendicular baseline.
%   HA = HEIGHTAMB(BPERP) returns height ambiguities of matrix BPERP.
%   Fixed parameters: wavelength=5.6 cm, slantrange=850 km, theta=21.5 deg.
%

% changed 18.5 deg to 21 deg for center pixel 09-Jan-2001 BK.

% $Revision: 1.2 $  $Date: 2001/03/16 13:46:42 $
% Bert Kampes 14-Dec-2000

if(nargin==0) help('heightamb'); error('no input'); end; 
lambda = 5.66e-2;
range1 = 850e3;
theta  = deg2rad(21.5);
ha     = (lambda.*range1.*sin(theta))./(2*bperp);

%%% EOF
