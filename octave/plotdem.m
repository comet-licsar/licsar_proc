% PLOTDEM - script to plot DEM from lat/lon/hei data.
%    Data is assumed to be in 3 files with names in the strings:
%    LATFILE, LONFILE, and HEIFILE (lowercase).
%    These variables should be in the workspace.
%    The file format is assumed to be real4 (as by Doris software).
%    You are prompted to give the number of lines in the files.
%    The files are read in matrices LAT, LON, HEI.
%    If these matrices exist in the workspace,
%    then they are not read in again.
%    Data equal to -999 is set to NaN (as by Doris software).
%
%    Basically a SURF plot is made, with handle H to it.
%
%    See also SURF, IMAGESC, FREADBK, VIEW

% $Revision: 1.3 $  $Date: 2000/06/28 11:52:45 $
% Bert Kampes, 28-Mar-2000



%%% Initialize.
more off;
disp('PLOTDEM: not tested...');
if (~(exist('LAT','var') & exist('LON','var') & exist('HEI','var')))
  if (~(exist('latfile','var') & exist('lonfile','var') & exist('heifile','var')))
    error('data not in workspace and no file names specified.');
  end;
  if (~exist(latfile,'file')) errror('latfile not found.'); end;
  if (~exist(lonfile,'file')) errror('lonfile not found.'); end;
  if (~exist(heifile,'file')) errror('heifile not found.'); end;
  % Read in data
  if (~exist('numlines','var'))
    numlines = input('Enter number of lines in the files: ');
  end;
  disp(['Reading matrix LAT from: ',latfile,'...']);
  LAT = freadbk(latfile,numlines,'float32'); 
  disp(['Reading matrix LON from: ',lonfile,'...']);
  LON = freadbk(lonfile,numlines,'float32'); 
  disp(['Reading matrix HEI from: ',heifile,'...']);
  HEI = freadbk(heifile,numlines,'float32'); 
  %%% remove artifical NaN's from matrices
  disp('Removing NaN''s from data...');
  xxx=find(HEI==-999);
  HEI(xxx)=NaN; LAT(xxx)=NaN; LON(xxx)=NaN;
  clear xxx;
end;


%%% Actually plot DEM.
disp('Plotting surface...');
h=surf(LON,LAT,HEI);
set(h,'edgecolor','none');
axis tight;
view(2);
colorbar;
title  'DEM'
xlabel 'geographic longitude (lambda)'
ylabel 'geographic latitude (phi)'

disp('saving in ps might be too large, rather jpg or grab it.');
disp('All done!');

% tidy up
more on;
%EOF
