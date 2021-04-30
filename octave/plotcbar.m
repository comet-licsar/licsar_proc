% PLOTCBAR -- Plot colorbar for current colormap.
%
%   PLOTCBAR is a script to plot a (vertical) colorbar for the
%   current colormap in figure 1. and saves it to a level 1 encapsulated
%   postscript file named cbar.eps
%   This eps can be used in latex documents.
%
%   To plot the DEOS2 colorbar, use:
%     colormap(deos2(256));
%     plotcbar;
%
%   See also PH, DEOS, DEOS2, JET, HOT, HSV
%

%// Bert Kampes, 26-Jan-2001

%%% create linear vector to imagesc.
np = 256;
a  = -pi:(2*pi)/(np-1):pi;
a = flipud(a(:));
figure(1);
h  = imagesc(a);

%%% Set size, ticks, ?
xmin  = 0.15;
ymin  = 0.15;
xsize = 0.07;
ysize = 0.7;
%set(1,'paper:...)
set(get(h,'parent'),'Position',[xmin ymin xsize ysize]);
set(get(h,'parent'),'xtick',[]);
set(get(h,'parent'),'xticklabel',[]);
set(get(h,'parent'),'ytick',[]);
set(get(h,'parent'),'yticklabel',[]);
%ticks=[-pi,0,pi];
%set(get(h,'parent'),'ytick',ticks);
%set(get(h,'parent'),'yticklabel',ticks);

%%% Add text.
t1=text(1,np*1.05,'-\pi');% center -pi on [0.5:1.5]
t2=text(1,-0.05*np,'\pi');% center  pi on [0.5:1.5]
set(t1,'fontname','helvetica');
set(t1,'fontunits','points');
set(t1,'fontsize',12);
set(t1,'fontangle','normal');
set(t1,'fontweight','normal');
set(t1,'HorizontalAlignment','center');
set(t2,'fontname','helvetica');
set(t2,'fontunits','points');
set(t2,'fontsize',12);
set(t2,'fontangle','normal');
set(t2,'fontweight','normal');
set(t2,'HorizontalAlignment','center');


%%% generate eps to use in latex
print -f1 -depsc cbar.eps

%%% EOF

