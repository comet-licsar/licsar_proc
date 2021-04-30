% InSAR utility toolbox.
% Version 1.3  26-Jan-2001
%
% MATLAB Toolbox for InSAR -- Delft University of Technology.
% 
%
% InSAR Utilities:
%   cpxmultilook    - Multilooking (spatial averaging) of complex data.
%   cpxdetrend      - Detrend complex interferogram.
%   heightamb       - Height ambiguity for given perpendicular baseline.
%   mirm            - Multi-Image Reflectifity Map.
%   myhamming       - Hamming window for filtering.
%   myrect          - Rectangular window for filtering.
%   oversample      - Oversampling of data by Fourier method (harmonic).
%   pdf_phase       - PDF for interferometric phase from coherence and multilook.
%   residues        - Locate residues in phase matrix.
%   std_phase       - Single look standard deviation based on coherence.
%   std_phasePS     - Phase standard deviation for point scatterers.
%   wl_phase_filt   - Wavelet phase filter, not included. Email request to:
%                     braunisc@mit.edu when interested.
%   wrap            - Wrap real values to principal interval [-pi,pi).
%
% Simulation:
%   bowl            - Generate unit bowl.
%   cone            - Generate unit cone.
%   pyramid         - Generate unit pyramid.
%   ramp            - Generate unit ramp.
%   siminterf       - Simulate phase of interferogram (nice!).
%   simnoise        - Simulate noise (for siminterf).
%   simstack        - Simulate datacube with phase (use approximations).
%   simslc          - Simulate SLC image.
%
% Unwrapping Tools:
%   manuwbk         - Manunally unwrap regions (tie, gui).
%   manuwbk_cb      - Call back functions for manuwbk.
%
% Data Utilities:
%   cc              - Clears workspace and closes all plots.
%   datatypesize    - Return number of bytes used by type.
%   freadbk         - Read a binary matrix on disk into matrix (nice!).
%   freadhgt        - Read a binary matrix on disk (hgt format).
%   freadras        - Read a SUN raster file from disk.
%   fsize           - Return number of bytes of file.
%   fwritebk        - Write a matrix to disk in binary format.
%   fwritehgt       - Write hgt matrix to disk (hgt format).
%   helphelp        - Called if input argument error to call help.
%   iseven          - Return 1 for even numbers, 0 otherwise.
%   isint           - Return 1 for integers, 0 otherwise.
%   isodd           - Return 1 for odd numbers, 0 otherwise.
%   ispow2          - Return 1 for powers of two, 0 otherwise.
%   issamesize      - Return 1 for arrays of same size, 0 otherwise.
%   isscalar        - Return 1 for scalars, 0 otherwise.
%   isvector        - Return 1 for vectors, 0 otherwise.
%   lying           - Returns lying vector.
%   normalize       - Rescale data vector to new interval.
%   standing        - Returns standing vector.
%
% Plot Utilities:
%   cpt2map         - Convert cpt table (GMT) to colormap.
%   deos            - Colormap for interferometric phase.
%   deos2           - Colormap for interferometric phase.
%   enlargefig      - Enlarge papersize of figure by 2 for printing.
%   mph             - Magnitude/phase display of complex interferogram.
%   ph              - Colormap for interferometric phase.
%   plotcbar        - Create eps for a single colorbar for phase.
%   plotdem         - Script to plot DEM from lat/lon/hei data files.
%   pra4            - Print figure fitted to A4 paper size.
%   titlebold       - Set figure title property to bold.
%   titlenormal     - Reset figure title.
%   tip             - Position figures over whole screen.
%   top             - Position figures on top of eachother.
%   trap            - Position figures slightly shifted wrt. eachother.
%
%
% See also FRACTAL, IMAGES, INSARDEMOS
% Tested on Matlab 5.3.0.10183 (R11) Jan 21 1999
% Tested on Matlab 6.1.0.450 (R12.1) May 18 2001
%

%%% EOF
%%% $Revision: 1.26 $  $Date: 2001/09/28 14:24:26 $
%%% Bert Kampes, 03-Mar-2000


