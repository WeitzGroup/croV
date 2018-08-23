%calculate burst size for each biological replicate, and group for table
T = getdatatable;

T_at15 = T(T.time==15,:);
T_at23 = T(T.time==23,:);

delV = T_at23.vfcm-T_at15.vfcm;
delH = T_at23.hosts-T_at15.hosts;
%std with propagated uncertainty
delHsd = sqrt(sum([T_at23.hosts_sd T_at15.hosts_sd].^2,2));
%burst size
burstsize = -delV./delH;
%std with propagated uncertainty -- no V uncertainty because no technical
%replicates of FCM
burstsizesd = sqrt(sum([(delV./(delH.^2).*delHsd)].^2,2));

%table -- group by MOI,
lowburst = mean(burstsize(4:5))
lowburstsd = std_biotech(lowburst,burstsizesd(4:5))

highburst = mean(burstsize(8:9))
highburstsd = std_biotech(highburst,burstsizesd(8:9))

[~,burstpvalue] = ttest2(burstsize(4:5),burstsize(8:9))