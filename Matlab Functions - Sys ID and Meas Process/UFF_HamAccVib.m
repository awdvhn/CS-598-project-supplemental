function [ time, F_Hammer, acc_Accelerometer, vel_Vibrometer ] = UFF_HamAccVib(  uff_Names )

uffRead_Cell = readuff(uff_Names);

time = uffRead_Cell{1}.x';
F_Hammer = uffRead_Cell{1}.measData';
acc_Accelerometer = detrend(uffRead_Cell{2}.measData)';
vel_Vibrometer = detrend(uffRead_Cell{3}.measData)';

end

