function hvalstruct = timing_pvalue(T,currdata,ystring,tlag)
%plots the one-sided p-values using paired t-tests for growth between data
%points separated by tlag (tlag refers to distance between points, which
%may or may not correspond to time depending on spread of data)

%MOI=1
[hvals,pvals,tvec] = subrun(T,currdata,1,tlag);
hvalstruct.MOI1 = hvals;
plot(tvec,pvals,'-o','LineWidth',2)
%MOI=10
[hvals,pvals,tvec] = subrun(T,currdata,10,tlag);
hvalstruct.MOI10 = hvals;
hold on
plot(tvec,pvals,'-o','LineWidth',2)
xl = xlim;
plot(xl,.05*ones(size(xl)),'k--','LineWidth',2)
hold off
xlabel('hours post-infection')
fullystring = strcat(ystring,{' growth p-value'});
ylabel(fullystring)
legend('MOI = 1','MOI = 10','Location','SouthEast')
set(gca,'Yscale','Log','FontSize',15)
end

function [hvals,pvals,tvec] = subrun(T,currdata,currMOI,tlag)
tearly = [0 1 2 3 4 5 6 7 8 9 12];
tlate = [12 13 14 15 16 17 18 19 20 21 22 23 24];
%set x value of plotted points at the average of the time for the data used
%in t-tests
tvec = ([tearly(1:end-tlag) tlate(1:end-tlag)]+...
    [tearly(tlag+1:end) tlate(tlag+1:end)])./2;
%pre-allocated significance vector
hvals = zeros(size(tvec));
%pre-allocate pvalue vector
pvals = zeros(size(tvec));
cnt = 0;
%loop through all points that have a corresponding point tlag away
%early experiment points
for jj = 1:(length(tearly)-tlag)
    cnt = cnt+1;
    beforelogic = T.time==tearly(jj) & T.MOI==currMOI & T.isearly==1;
    afterlogic = T.time==tearly(jj+tlag) & T.MOI==currMOI & T.isearly==1;
    beforevals = currdata(beforelogic);
    aftervals = currdata(afterlogic);
    %one sided paired t-test
    [hvals(cnt),pvals(cnt)] = ttest(beforevals,aftervals,'Alpha',.05,'Tail','Left');
end

%loop through all points that have a corresponding point tlag away
%late experiment points
for jj = 1:(length(tlate)-tlag)
    cnt = cnt+1;
    beforelogic = T.time==tlate(jj) & T.MOI==currMOI & T.isearly==0;
    afterlogic = T.time==tlate(jj+tlag) & T.MOI==currMOI & T.isearly==0;
    beforevals = currdata(beforelogic);
    aftervals = currdata(afterlogic);
    %one sided paired t-test
    [hvals(cnt),pvals(cnt)] = ttest(beforevals,aftervals,'Alpha',.05,'Tail','Left');
end
end