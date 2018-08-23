function plot2pdf(filename,varargin)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
if nargin == 1
ax_width = outerpos(3) - ti(1) - ti(3);
else
    ax_width = .975*(outerpos(3) - ti(1) - ti(3));
end
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
strname = strcat('figures/',filename);
print(fig,strname,'-dpdf')