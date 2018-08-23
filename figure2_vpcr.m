function [vinit,vinf,vfin] = figure2_vpcr(MOI,varargin)
close all
fsize=20;
msize=30;
ds = getdata;%structure with data
moi1.early.rep1 = ds.vpcr(ds.MOI==MOI & ds.isearly==1 & ds.replicate==1);
moi1.early.rep2 = ds.vpcr(ds.MOI==MOI & ds.isearly==1 & ds.replicate==2);
moi1.early.rep3 = ds.vpcr(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);
moi1.early.time = ds.time(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);
moi1.early.rep1sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==1);
moi1.early.rep2sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==2);
moi1.early.rep3sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);

moi1.late.rep1 = ds.vpcr(ds.MOI==MOI & ds.isearly==0 & ds.replicate==1);
moi1.late.rep2 = ds.vpcr(ds.MOI==MOI & ds.isearly==0 & ds.replicate==2);
moi1.late.rep3 = ds.vpcr(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);
moi1.late.time = ds.time(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);
moi1.late.rep1sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==1);
moi1.late.rep2sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==2);
moi1.late.rep3sd = ds.vpcr_sd(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);

vinit = [moi1.early.rep1(1) moi1.early.rep2(1) moi1.early.rep3(1) ...
    moi1.late.rep1(1) moi1.late.rep2(1) moi1.late.rep3(1)];
vinf = [moi1.early.rep1(2) moi1.early.rep2(2) moi1.early.rep3(2) ...
    moi1.late.rep1(2) moi1.late.rep2(2) moi1.late.rep3(2)];
vfin = [moi1.early.rep1(end) moi1.early.rep2(end) moi1.early.rep3(end) ...
    moi1.late.rep1(end) moi1.late.rep2(end) moi1.late.rep3(end)];

figure(1)
plot(moi1.early.time,moi1.early.rep1,'x','MarkerSize',msize-18,'LineWidth',2)
hold on
plot(moi1.early.time,moi1.early.rep2,'x','MarkerSize',msize-18,'LineWidth',2)
plot(moi1.early.time,moi1.early.rep3,'x','MarkerSize',msize-18,'LineWidth',2)

plot(moi1.late.time,moi1.late.rep1,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep2,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep3,'.','MarkerSize',msize)

ax = gca;
ax.ColorOrderIndex = 1;
errorbar(moi1.early.time,moi1.early.rep1,moi1.early.rep1sd,'x','MarkerSize',msize-18,'LineWidth',2)
errorbar(moi1.early.time,moi1.early.rep2,moi1.early.rep2sd,'x','MarkerSize',msize-18,'LineWidth',2)
errorbar(moi1.early.time,moi1.early.rep3,moi1.early.rep3sd,'x','MarkerSize',msize-18,'LineWidth',2)

errorbar(moi1.late.time,moi1.late.rep1,moi1.late.rep1sd,'.','MarkerSize',msize,'LineWidth',2)
errorbar(moi1.late.time,moi1.late.rep2,moi1.late.rep2sd,'.','MarkerSize',msize,'LineWidth',2)
errorbar(moi1.late.time,moi1.late.rep3,moi1.late.rep3sd,'.','MarkerSize',msize,'LineWidth',2)
hold off
set(gca,'FontSize',fsize)
ylim([0 15*10^8])
if MOI == 1
title('low MOI (1)','FontSize',fsize)
else
    title('high MOI (10)','FontSize',fsize)
end
xlabel('hpi','FontSize',fsize)
ylabel('Viral DNA copies per ml','FontSize',fsize)
xlim([0 25])
%if nargin==2
%    [h,objs] = columnlegend(2,{'E1','E2','E3','L1','L2','L3'},'Location','SouthEast','FontSize',fsize);
%    set(findall(objs, 'type', 'text'), 'fontsize', fsize)
%    set(gca,'YScale',varargin{1},'FontSize',fsize)
%else
%    legend('E1','E2','E3','L1','L2','L3','Location','NorthWest')
%set(gca,'YScale','Linear','FontSize',fsize)
%end

[h,objs] = columnlegend(2,{'E1','E2','E3','L1','L2','L3'},'Location','SouthEast','FontSize',fsize);
set(findall(objs, 'type', 'text'), 'fontsize', fsize)
set(gca,'YScale',varargin{1},'FontSize',fsize)

try
    plot2pdf(strcat('vpcr_',num2str(MOI),varargin{1}))
catch
    plot2pdf(strcat('vpcr_',num2str(MOI)))
end