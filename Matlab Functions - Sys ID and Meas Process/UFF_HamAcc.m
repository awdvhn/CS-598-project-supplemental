function [ time, F_Hammer, acc_Accelerometer1, acc_Accelerometer2] = UFF_HamAcc(  uff_File )

uffRead_Cell = readuff(uff_File);

time = uffRead_Cell{1}.x';
acc_Accelerometer1 = detrend(uffRead_Cell{1}.measData)';
acc_Accelerometer2 = detrend(uffRead_Cell{2}.measData)';

F_Hammer = uffRead_Cell{3}.measData';

end

