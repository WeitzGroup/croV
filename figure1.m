%set seed for consistency
rng(1337)   %h4x0rZ

%get data
dstable = getdatatable;

%Look at MOI=0, no viruses
dstable_MOI0 = dstable(dstable.MOI==0,:);   %only get control MOI=0
dstable_MOI0 = sortrows(dstable_MOI0,{'isearly','replicate','time'}); %order to eventually grab before and after 
tmin = 4;   %ignore times before clear batch effect loss
logicbefore = (dstable_MOI0.time>=tmin & dstable_MOI0.time<=9) & dstable_MOI0.isearly==1 | ...
    (dstable_MOI0.time>=12 & dstable_MOI0.time<=23) & dstable_MOI0.isearly==0;
logicafter = [false; logicbefore(1:end-1)];
%make ordered tables of before and after times
tblbefore = dstable_MOI0(logicbefore,:);
tblafter = dstable_MOI0(logicafter,:);

close all
%plot the raw data with  error bars
figure(1)
%font and marker size
fsize = 20;
msize = 30;

%set order of data to plot
dstable_MOI0.plotidx = dstable_MOI0.replicate+3*(1-dstable_MOI0.isearly);

%plot data in order of replicate to use consistent color and marking
uniplotidx = unique(dstable_MOI0.plotidx);
for jj = 1:length(uniplotidx)
    currlogic = dstable_MOI0.plotidx==uniplotidx(jj);
    if jj<4
        plot(dstable_MOI0.time(currlogic),dstable_MOI0.hosts(currlogic),'x','MarkerSize',msize-18,'LineWidth',2)
    else
        plot(dstable_MOI0.time(currlogic),dstable_MOI0.hosts(currlogic),'.','MarkerSize',msize)
    end
    hold on
end
legend('CE1','CE2','CE3','CL1','CL2','CL3','Location','NorthWest','AutoUpdate','off');
legend boxoff
%plot with error bar (first loop was too set the legend, cos matlab can be
%annoying)
ax = gca;
ax.ColorOrderIndex = 1;
for jj = 1:length(uniplotidx)
    currlogic = dstable_MOI0.plotidx==uniplotidx(jj);
    if jj<4
    errorbar(dstable_MOI0.time(currlogic),dstable_MOI0.hosts(currlogic),dstable_MOI0.hosts_sd(currlogic),...
    'x','MarkerSize',msize-18,'LineWidth',2)
    else
        errorbar(dstable_MOI0.time(currlogic),dstable_MOI0.hosts(currlogic),dstable_MOI0.hosts_sd(currlogic),...
    '.','MarkerSize',msize,'LineWidth',2)
    end
    hold on
end
hold off
xlim([-1 25])
xlabel('Time (hours)','FontSize',fsize)
ylabel('Host cells per ml','FontSize',fsize)
set(gca,'YScale','Linear','FontSize',fsize)
plot2pdf('control_lin')


%Now plot regression on averaged growth rates.

%The host abundances are manually counted and can fluctuate due to
%subsampling, which is exacerbated by spatial heterogeneity in the flask.
%We average host abundances across biological replicates for each time
%point to estimate growth rate.

%group
[unisamp,ia] = unique(tblbefore.sample,'stable');
%3 cos of number of replicates
grptblbefore = tblbefore(1:3:(1+3*(length(unisamp)-1)),:);
pairtbl = tblbefore(ia,:);
%preallocate new rows
pairtbl.hosts_after = zeros(size(pairtbl.hosts));
pairtbl.hosts_after_sd = zeros(size(pairtbl.hosts_sd));
pairtbl.tau = zeros(size(pairtbl.time));
for jj= 1:length(unisamp)
    beforesamp = unisamp(jj);
    aftersamp = beforesamp+3;%sampling convention, see data file
    currbeforelogic = dstable_MOI0.sample==beforesamp;
    pairtbl.hosts(jj) = mean(dstable_MOI0.hosts(currbeforelogic));
    %for stdev propagate uncertainty using std_biotech
    pairtbl.hosts_sd(jj) = std_biotech(dstable_MOI0.hosts(currbeforelogic),dstable_MOI0.hosts_sd(currbeforelogic));
    timebefore = mean(dstable_MOI0.time(currbeforelogic)); %mean just to turn vector into scalar
    currafterlogic = dstable_MOI0.sample==aftersamp;
    pairtbl.hosts_after(jj) = mean(dstable_MOI0.hosts(currafterlogic));
    %for stdev propagate uncertainty using std_biotech
    pairtbl.hosts_after_sd(jj) = std_biotech(dstable_MOI0.hosts(currafterlogic),dstable_MOI0.hosts_sd(currafterlogic));
    timeafter = mean(dstable_MOI0.time(currafterlogic)); %mean just to turn vector into scalar
    pairtbl.tau(jj) = timeafter-timebefore;
end


currx = pairtbl.hosts;
currxsd = pairtbl.hosts_sd;
%growth rate: y= d(log(x))/dt, use numerical approximation
curry = (log(pairtbl.hosts_after)-log(pairtbl.hosts))./pairtbl.tau;
%propagate uncertainty for calculating standard deviation
currysd = sqrt(sum((1./pairtbl.tau.*[1./pairtbl.hosts_after.*pairtbl.hosts_after_sd ...
    -1./pairtbl.hosts.*pairtbl.hosts_sd]).^2,2));

figure(2)
%this function fits line to the data accounting for the stdev of each
%data point
 [p,psd] = linfitxy(currx,curry,currxsd,currysd...
     ,'Plotting',false);
 plot(currx,curry,'b.','MarkerSize',msize)
 
 currfun = @(x) p(2)+p(1)*x;
 xrange = [min(currx-currxsd) max(currx+currxsd)];
 hold on
 plot(xrange,currfun(xrange),'r','LineWidth',1.5)
 errorbar(currx,curry,currysd,currysd,currxsd,currxsd,'b.','MarkerSize',msize)
 hold off
 xlim(xrange)
 legend('Data','Linear Fit')
 title('')
xlabel('Host (H) abundance')
ylabel('$\frac{\Delta \log H}{\Delta t}$','interpreter','Latex')
set(gca,'FontSize',18,'Xscale','Linear')
plot2pdf('control_linfit',1)

%pars: [r K], growth rate and carrying capacity
pars = [p(2) p(2)./p(1)];
rK = pars
%propagate uncertainty for stdev
pars_sd = [psd(2) sqrt(sum([(-p(2)/(p(1)^2))*psd(1) (1/p(1))*psd(2)].^2))];
%calculate doubling time
doublingtime = log(2)./pars(1)
%approximate stdev of a function using Taylor expansion: 
%Var(f(x)) = Var(x)*(f'(X))^2
doublingtime_sd = sqrt((pars_sd(1).^2)*((-log(2)./(pars(1).^2)).^2))