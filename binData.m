function [mean_bin, centers] = binData(activity, behaviour, edges)

% returns for each bin, specificed by the edges which bin each data index
% lies in (bin) as well as the number of data points in each bin (N)
[N, edges, bin] = histcounts(behaviour(:,1), edges);
% sums datapoints belonging to the same bin & saves the values in temp
temp = accumarray(bin+1, activity, [length(edges) 1]);
% finds the average of each bin value
mean_bin = bsxfun(@rdivide, temp(2:end), N');
% finds the center value of each bin
centers = edges(1:end-1)+diff(edges)/2;

end