function figure5(numboxcar)

%function to plot figure 5 simulated dynamics. numboxcar sets the number of
%infection subcompartments as shown in title of subfigures.
close all
T = getdatatable;
%min latency time to start simulation
lattime = 6;


currMOI=1;
%initial conditions based on data
Htot = mean(T.hosts(T.time==6 & T.MOI==currMOI));
%rough estimate number of infected based on poisson distribution
Hinf = .5*mean(T.hosts(T.time==2 & T.MOI==currMOI));
[sol1] = myddemodel(lattime,Htot,Hinf,numboxcar);

currMOI=10;
%initial conditions based on data
Htot = mean(T.hosts(T.time==6 & T.MOI==currMOI));
%rough estimate number of infected based on poisson distribution
Hinf = .9*mean(T.hosts(T.time==2 & T.MOI==currMOI));
[sol10] = myddemodel(lattime,Htot,Hinf,numboxcar);

figure(1)
plot(sol1.x+lattime,sol1.y(1,:),'LineWidth',2)
hold on
plot(sol10.x+lattime,sol10.y(1,:),'LineWidth',2)
yl = [0 4*10^8];
plot(lattime*ones(size(yl)),yl,'k--','LineWidth',2)
hold off
xlabel('hours post-infection')
ylabel('Viruses per ml')
ylim(yl)
xlim([0 25])
legend('MOI = 1','MOI=10','Location','NorthWest')
set(gca,'FontSize',18)
filestr = strcat('modelv_n',num2str(numboxcar));
plot2pdf(filestr)

figure(2)
plot(sol1.x+lattime,sum(sol1.y(2:end,:),1),'LineWidth',2)
hold on
plot(sol10.x+lattime,sum(sol10.y(2:end,:),1),'LineWidth',2)
yl = [0 1*10^6];
plot(lattime*ones(size(yl)),yl,'k--','LineWidth',2)
hold off
xlabel('hours post-infection')
ylabel('Hosts per ml')
%title('multiple stages of infection')
if numboxcar==1
    title(['n = ' num2str(numboxcar) ' stage of infection'])
else
title(['n = ' num2str(numboxcar) ' stages of infection'])
end
ylim(yl)
xlim([0 25])
legend('MOI = 1','MOI=10')
set(gca,'FontSize',18)
filestr = strcat('modelh_n',num2str(numboxcar));
plot2pdf(filestr)