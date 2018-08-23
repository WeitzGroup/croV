function [inithosts,minhosts,delh2,predelh2] = figure2_hosts(MOI,varargin)
close all
fsize=20;
msize = 30;
ds = getdata;%structure with data
moi1.early.rep1 = ds.hosts(ds.MOI==MOI & ds.isearly==1 & ds.replicate==1);
moi1.early.rep2 = ds.hosts(ds.MOI==MOI & ds.isearly==1 & ds.replicate==2);
moi1.early.rep3 = ds.hosts(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);
moi1.early.time = ds.time(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);
moi1.early.rep1sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==1);
moi1.early.rep2sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==2);
moi1.early.rep3sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);

moi1.late.rep1 = ds.hosts(ds.MOI==MOI & ds.isearly==0 & ds.replicate==1);
moi1.late.rep2 = ds.hosts(ds.MOI==MOI & ds.isearly==0 & ds.replicate==2);
moi1.late.rep3 = ds.hosts(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);
moi1.late.time = ds.time(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);
moi1.late.rep1sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==1);
moi1.late.rep2sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==2);
moi1.late.rep3sd = ds.hosts_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);

control.early.rep1 = ds.hosts(ds.MOI==0 & ds.isearly==1 & ds.replicate==1);
control.early.rep2 = ds.hosts(ds.MOI==0 & ds.isearly==1 & ds.replicate==2);
control.early.rep3 = ds.hosts(ds.MOI==0 & ds.isearly==1 & ds.replicate==3);
control.early.time = ds.time(ds.MOI==0 & ds.isearly==1 & ds.replicate==3);

meancontrol.early = (control.early.rep1+control.early.rep2+control.early.rep3)/3;

control.late.rep1 = ds.hosts(ds.MOI==0 & ds.isearly==0 & ds.replicate==1);
control.late.rep2 = ds.hosts(ds.MOI==0 & ds.isearly==0 & ds.replicate==2);
control.late.rep3 = ds.hosts(ds.MOI==0 & ds.isearly==0 & ds.replicate==3);
control.late.time = ds.time(ds.MOI==0 & ds.isearly==0 & ds.replicate==3);

meancontrol.late = (control.late.rep1+control.late.rep2+control.late.rep3)/3;

inithosts = [moi1.early.rep1(1) moi1.early.rep2(1) moi1.early.rep3(1) ...
    moi1.late.rep1(1) moi1.late.rep2(1) moi1.late.rep3(1)];
minhosts = [min(moi1.early.rep1(1:4)) min(moi1.early.rep2(1:4)) min(moi1.early.rep3(1:4)) ...
    min(moi1.late.rep1(1:2)) min(moi1.late.rep2(1:2)) min(moi1.late.rep3(1:2))];
delh2 = [moi1.late.rep1(moi1.late.time==14)-moi1.late.rep1(moi1.early.time==24) ...
    moi1.late.rep2(moi1.late.time==14)-moi1.late.rep2(moi1.early.time==24) ...
    moi1.late.rep3(moi1.late.time==14)-moi1.late.rep3(moi1.early.time==24)];

predelh2 = [moi1.late.rep1(moi1.late.time==14) ...
    moi1.late.rep2(moi1.late.time==14) ...
    moi1.late.rep3(moi1.late.time==14)];
figure(1)
plot(moi1.early.time,moi1.early.rep1,'x','MarkerSize',msize-18,'LineWidth',2)
hold on
plot(moi1.early.time,moi1.early.rep2,'x','MarkerSize',msize-18,'LineWidth',2)
plot(moi1.early.time,moi1.early.rep3,'x','MarkerSize',msize-18,'LineWidth',2)

plot(moi1.late.time,moi1.late.rep1,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep2,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep3,'.','MarkerSize',msize)
plot(moi1.early.time,meancontrol.early,'kx','MarkerSize',msize-18,'LineWidth',2)
plot(moi1.late.time,meancontrol.late,'.','Color',[169 169 169]./255,'MarkerSize',msize)

ax = gca;
ax.ColorOrderIndex = 1;
errorbar(moi1.early.time,moi1.early.rep1,moi1.early.rep1sd,'x','MarkerSize',msize-18,'LineWidth',2)
errorbar(moi1.early.time,moi1.early.rep2,moi1.early.rep2sd,'x','MarkerSize',msize-18,'LineWidth',2)
errorbar(moi1.early.time,moi1.early.rep3,moi1.early.rep3sd,'x','MarkerSize',msize-18,'LineWidth',2)

errorbar(moi1.late.time,moi1.late.rep1,moi1.late.rep1sd,'.','MarkerSize',msize,'LineWidth',2)
errorbar(moi1.late.time,moi1.late.rep2,moi1.late.rep2sd,'.','MarkerSize',msize,'LineWidth',2)
errorbar(moi1.late.time,moi1.late.rep3,moi1.late.rep3sd,'.','MarkerSize',msize,'LineWidth',2)

hold off
legend('E1','E2','E3','L1','L2','L3','\langleEC\rangle','\langleLC\rangle','Location','NorthWest')
xlim([0 25])
if MOI == 1
title('low MOI (1)','FontSize',fsize)
else
    title('high MOI (10)','FontSize',fsize)
end
xlabel('hpi','FontSize',fsize)
ylabel('Hosts per ml','FontSize',fsize)
if nargin==2
    set(gca,'YScale',varargin{1},'FontSize',fsize)
    [h,objs] = columnlegend(3,{'E1','E2','E3','L1','L2','L3','\langleEC\rangle','\langleLC\rangle'}...
        ,'Location','SouthWest','FontSize',fsize);
    set(findall(objs, 'type', 'text'), 'fontsize', fsize)
else
    legend('E1','E2','E3','L1','L2','L3','\langleEC\rangle','\langleLC\rangle','Location','NorthWest')
set(gca,'YScale','Linear','FontSize',fsize)
end
plot2pdf(strcat('hosts_',num2str(MOI)))
