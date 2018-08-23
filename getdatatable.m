function T = getdatatable
%read in excel data and manually add column names

ds = getdata;

T = struct2table(ds);