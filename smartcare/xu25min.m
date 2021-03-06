function y = xu25min(x)

% xu25min - returns the min of the data points excluding the upper 25% in x (sorted by
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

percentile75 = round(size(x, dim) * .75);

y = min(x(1:percentile75));

end