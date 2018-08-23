T = getdatatable;
%focus on MOI=10 first
currMOIlogic = T.MOI==10;

%Calculate adsorption rate, phi, for y-values
V0 = T.vpcr(currMOIlogic & T.time==0);
V1 = T.vpcr(currMOIlogic & T.time==1);
V0sd = T.vpcr_sd(currMOIlogic & T.time==0);
V1sd = T.vpcr_sd(currMOIlogic & T.time==1);

H0 = T.hosts(currMOIlogic & T.time==0);
H1 = T.hosts(currMOIlogic & T.time==1);
H0sd = T.hosts_sd(currMOIlogic & T.time==0);
H1sd = T.hosts_sd(currMOIlogic & T.time==1);

[phi10,phisd10] = calculatephi(V0,V1,H0,H1,V0sd,V1sd,H0sd,H1sd);
ally = [];
allysd = [];
ally = [ally; phi10];
allysd = [allysd; phisd10];
meanphi10 = mean(phi10);
meanphi10sd = std_biotech(phi10,phisd10);

%Calculate realized MOI of innoculant
[currvh,currvhsd] =calculateMOI(V0,H0,V0sd,H0sd);
allx = [];
allxsd = [];
allx = [allx; currvh];
allxsd = [allxsd; currvhsd];

%focus on MOI=1 now
currMOIlogic = T.MOI==1;

V0 = T.vpcr(currMOIlogic & T.time==0);
V1 = T.vpcr(currMOIlogic & T.time==1);
V0sd = T.vpcr_sd(currMOIlogic & T.time==0);
V1sd = T.vpcr_sd(currMOIlogic & T.time==1);

H0 = T.hosts(currMOIlogic & T.time==0);
H1 = T.hosts(currMOIlogic & T.time==1);
H0sd = T.hosts_sd(currMOIlogic & T.time==0);
H1sd = T.hosts_sd(currMOIlogic & T.time==1);

[phi1,phisd1] = calculatephi(V0,V1,H0,H1,V0sd,V1sd,H0sd,H1sd);
ally = [ally; phi1];
allysd = [allysd; phisd1];
meanphi1 = mean(phi1);
meanphi1sd = std_biotech(phi1,phisd1);

[currvh,currvhsd] =calculateMOI(V0,H0,V0sd,H0sd);
allx = [allx; currvh];
allxsd = [allxsd; currvhsd];
close all
%handaxes1 = axes('Position', [0.12 0.12 0.8 0.8]); 
msize=30;
% Main plot -- log-log to see data better
%plot once to get axis limits
plot(allx,ally,'b.','MarkerSize',msize)
maxy = max(ally+allysd);
set(gca,'Yscale','Log','Xscale','Log','FontSize',16)
xl = xlim;
yl = ylim;
ylim([yl(1) 1.1*maxy])
yl = ylim;
%fit exponential to the data
mdl = fitlm(allx,log(ally));
%make function for plotting
mdlfun_beforetrans = @(x) mdl.Coefficients.Estimate(2).*x+mdl.Coefficients.Estimate(1);
%transform it back to original axes
mdlfun = @(x) exp(mdlfun_beforetrans(x));
%make domain points for plotting
xldomain = logspace(log(xl(1)),log(xl(2)),1000);
hold on
%plot fit
plot(xldomain,mdlfun(xldomain),'r','LineWidth',2)
%plot data with errorbars
errorbar(allx,ally,allysd,allysd,allxsd,allxsd,'b.','MarkerSize',msize)
hold off
legend('Data','Exponential Fit')
xlim(xl)
ylim(yl)
xlabel('Initial MOI')
ylabel('Adsorption rate, \phi')
set(gca,'Yscale','Log','Xscale','Log','FontSize',16)

% Adjust XY label font
handxlabel1 = get(gca, 'XLabel');
set(handxlabel1, 'FontSize', 16, 'FontWeight', 'bold')
handylabel1 = get(gca, 'ylabel');
set(handylabel1, 'FontSize', 16, 'FontWeight', 'bold')

% Inset -- semilog plot, same as above but different scaling on axes
handaxes2 = axes('Position', [0.2 0.2 0.2 0.2]);
%plot to get axes
plot(allx,ally,'b.','MarkerSize',msize-10)
maxy = max(ally+allysd);
set(gca,'Yscale','Log','Xscale','Linear','FontSize',6)
xl = xlim;
yl = ylim;
ylim([yl(1) 1.1*maxy])
yl = ylim;
xldomain = linspace(xl(1),xl(2),1000);
hold on
plot(xldomain,mdlfun(xldomain),'r','LineWidth',1)
errorbar(allx,ally,allysd,allysd,allxsd,allxsd,'b.','MarkerSize',msize-10)
hold off
xlim(xl)
ylim(yl)

set(handaxes2, 'Box', 'off')

% Adjust XY label font
set(get(handaxes2, 'XLabel'), 'FontName', 'Times')
set(get(handaxes2, 'YLabel'), 'FontName', 'Times')

plot2pdf('adsorption_scaling_inset')

%table phi
meanphi1
meanphi1sd
meanphi10
meanphi10sd
[~,phipvalue] = ttest2(phi1,phi10)