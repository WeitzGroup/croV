function [vinit,delv2] = figure2_vfcm(MOI,varargin)
close all
msize=30;
fsize = 20;
ds = getdata;%structure with data
moi1.early.rep1 = ds.vfcm(ds.MOI==MOI & ds.isearly==1 & ds.replicate==1);
moi1.early.rep2 = ds.vfcm(ds.MOI==MOI & ds.isearly==1 & ds.replicate==2);
moi1.early.rep3 = ds.vfcm(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);
moi1.early.time = ds.time(ds.MOI==MOI & ds.isearly==1 & ds.replicate==3);

moi1.late.rep1 = ds.vfcm(ds.MOI==MOI & ds.isearly==0 & ds.replicate==1);
moi1.late.rep2 = ds.vfcm(ds.MOI==MOI & ds.isearly==0 & ds.replicate==2);
moi1.late.rep3 = ds.vfcm(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);
moi1.late.time = ds.time(ds.MOI==MOI & ds.isearly==0 & ds.replicate==3);

delv2 = -[moi1.late.rep1(moi1.late.time==14)-max(moi1.late.rep1(moi1.early.time>12)) ...
    moi1.late.rep2(moi1.late.time==14)-max(moi1.late.rep2(moi1.early.time>12)) ...
    moi1.late.rep3(moi1.late.time==14)-max(moi1.late.rep3(moi1.early.time>12))];

vinit = [moi1.early.rep1(1) moi1.early.rep2(1) moi1.early.rep3(1) ...
    moi1.late.rep1(1) moi1.late.rep2(1) moi1.late.rep3(1)];
figure(1)
plot(moi1.early.time,moi1.early.rep1,'x','MarkerSize',msize-18,'LineWidth',2)
hold on
plot(moi1.early.time,moi1.early.rep2,'x','MarkerSize',msize-18,'LineWidth',2)
plot(moi1.early.time,moi1.early.rep3,'x','MarkerSize',msize-18,'LineWidth',2)
plot(moi1.late.time,moi1.late.rep1,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep2,'.','MarkerSize',msize)
plot(moi1.late.time,moi1.late.rep3,'.','MarkerSize',msize)
set(gca,'FontSize',fsize)
hold off

xlim([0 25])
ylim([0 7*10^8])
if MOI == 1
title('low MOI (1)','FontSize',fsize)
else
    title('high MOI (10)','FontSize',fsize)
end
xlabel('hpi','FontSize',fsize)
ylabel('Virus particles per ml','FontSize',fsize)
if nargin==2
    set(gca,'YScale',varargin{1},'FontSize',fsize)
else
    set(gca,'YScale','Linear','FontSize',fsize)
end

[h,objs] = columnlegend(2,{'E1','E2','E3','L1','L2','L3'},'Location','SouthEast','FontSize',fsize);
set(findall(objs, 'type', 'text'), 'fontsize', fsize)
set(gca,'YScale',varargin{1},'FontSize',fsize)
plot2pdf(strcat('vfcm_',num2str(MOI)))