close all
hvalstruct = timing_pvalue(T,T.vfcm,'CroV FCM',1);
plot2pdf('pval_fcm')
close all
hvalstruct = timing_pvalue(T,T.vpcr,'CroV qPCR',1);
plot2pdf('pval_pcr')