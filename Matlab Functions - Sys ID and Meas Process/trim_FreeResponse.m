function [t_trim, x_trim, trimRange] = trim_FreeResponse(tspan, x, x_minProp, cycles_to_skip)

%% First Trim: From max point to points where amplitude less than minProp of max 

% Find Highest Maximum 
[xMax, i_start1] = max(x);

% Find End where oscillations are smaller than min Proportion 
i_end = find(flip(x) > x_minProp * xMax, 5);
i_end = length(x) - i_end(end);

% First Trim
t_trim = tspan(i_start1:i_end) - tspan(i_start1);
x_trim = x(i_start1:i_end);

%% Second Trim: From maximum after skipped cycles to last maximum

% Skip the initial cycles required to be skipped
[f_trim, Xamp_trim] = regFFT(t_trim, x_trim);
[~, iMax] = max(Xamp_trim);
[~, i_start2] = max(x_trim(t_trim > (cycles_to_skip-0.5)/f_trim(iMax)));
i_start2 = i_start2 + find(t_trim > (cycles_to_skip-0.5)/f_trim(iMax), 1, 'First') - 1;

% End Trimming at the maximum
[~, i_end] = max(x_trim(t_trim > t_trim(end) - 1.25/f_trim(iMax)));
i_end = i_end + find(t_trim > t_trim(end) - 1.25/f_trim(iMax), 1, 'First') - 1;

% Final Trim
t_trim = t_trim(i_start2:i_end) - tspan(i_start2);
x_trim = x_trim(i_start2:i_end);
i_start = i_start1+i_start2-1;
i_end = i_end + i_start1-1;
trimRange = i_start:i_end;
end

