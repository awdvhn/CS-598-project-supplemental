function [F_DeNoised, i_Fstart, i_Fend] = hammerForce_trimDeNoise(F_raw)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

frac_Fstart = 0.5e-2; % ratio of maximum force to be considered as zero
N = length(F_raw);
F_DeNoised = [F_raw(1:floor(N/50));...
    F_raw(floor(N/50))*ones(N - floor(N/50),1)];
F_DeNoised = detrend(F_DeNoised, 0);

[Fmax, i_Fmax] = max(F_DeNoised);
i_Fstart = i_Fmax + 1 - find(F_DeNoised(i_Fmax:-1:1) <= frac_Fstart*Fmax, 1, 'First');
if F_DeNoised(i_Fstart) < 0
    F_DeNoised(i_Fstart) = 0;
end

i_Fend = i_Fmax - 1 + find(F_DeNoised(i_Fmax:end) <= frac_Fstart*Fmax, 1, 'First');
if F_DeNoised(i_Fend) < 0
    F_DeNoised(i_Fend) = 0;
end

F_DeNoised([1:i_Fstart-1 i_Fend+1:end]) = 0;

end

