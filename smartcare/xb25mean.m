function y = xb25mean(x)

% xb25mean - returns the mean of the data points excluding the bottom 25% in x (sorted by
% magnitude)  

x = sort(x, 'ascend');

if size(x, 1) == 1
    dim = 2;
elseif size(x, 2) == 1
    dim = 1;
else
    fprintf('**** function only works for vectors ****');
    return
end

percentile25 = round(size(x, dim) * .25) + 1;

y = mean(x(percentile25:end));

end

