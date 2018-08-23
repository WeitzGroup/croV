function datastruct = getdata
%read in excel data and manually add column names

filename = 'data/alldata_clean.xlsx';
[num,~,~] = xlsread(filename);
datastruct.sample = num(:,1);
datastruct.replicate = num(:,2);
datastruct.hosts = num(:,3);
datastruct.time = num(:,4);
datastruct.MOI = num(:,5);
datastruct.vfcm = num(:,6);
datastruct.vpcr = num(:,7);
datastruct.isearly = num(:,8);
datastruct.hosts_sd = num(:,9);
datastruct.vpcr_sd = num(:,10);