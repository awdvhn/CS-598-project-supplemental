function [ Cres_acc, R2res_acc, Pres_acc, tRes_acc, interSec0, interSec1 , interSec2] = ...
    compute_InitialCoeffGuess( m, t0, v0, time, data_2_zero, data1, data2,...
    var2_FUN_var1, C0, C_lb, C_ub)

% Get Points for Restoring Effect: xdot = 0
[tRes_acc, interSec0]= intersections(t0, v0,...
    time, data_2_zero); % Find instants of intersection b/w data0 and v0
interSec1 = interp1(time, data1, tRes_acc);% data1 of points at intersection points
interSec2 = interp1(time, data2, tRes_acc);% data2 of points at intersection points

% Fit Restoring points to compute initial Guesses
options = psoptimset;
options.Display = 'iter';
options.TolMesh = 1e-10;
options.MaxIter = 1e10;

objFun = @(C) sum((interSec2-var2_FUN_var1(C, m, interSec1)).^2)/sum((interSec2-mean(interSec2)).^2);

[Cres_acc, Pres_acc] = patternsearch(objFun,...
    C0,[],[],[],[],C_lb,C_ub,[],options);
R2res_acc = rsquare(interSec2, var2_FUN_var1(Cres_acc, m, interSec1));

end

