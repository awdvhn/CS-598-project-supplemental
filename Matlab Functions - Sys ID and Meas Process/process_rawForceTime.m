function [time, F_deNoise, iFstart, iFend,...
    timeForce] = process_rawForceTime(timeRAW, Fraw_Hammmer)

%%% Denoise Hammer Force, Get Impulse Start and End instants after denoise
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(Fraw_Hammmer); 

%%% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

%%% Trim Time Instants to after Hit Region
timeForce = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

%%% Trim Denoised Force to after Hit Region
F_deNoise = F_deNoise(I_ForceStart);

%% Trim Time
I_ForceEnds = iFend:length(timeRAW);
time = timeRAW(I_ForceEnds) - timeRAW(I_ForceEnds(1));

end

