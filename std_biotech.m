function stds_prop = std_biotech(mus,stds)
%standard deviation between datapoints that individually have std.
%refer to: https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
vars = stds.^2;
stds_prop = sqrt(mean(vars)+mean(mus.^2) - (mean(mus).^2));