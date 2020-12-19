function dX = DynSys_DoubleLatticesLZT(t, X, M, K, D, alpha_ALL,...
    beta_ALL, t_exp, Fexp_cell, t_Fend)
DoF_conf = size(M, 1);
%% Forcing to be Applied
Fapp_t = zeros(DoF_conf, 1);
if t < t_Fend
    for i_Forcing = 1:length(Fexp_cell)
        if isempty(Fexp_cell{i_Forcing}) ~= 1
            Fapp_t(i_Forcing) = interp1(t_exp, Fexp_cell{i_Forcing}, t);
        else
            Fapp_t(i_Forcing) = 0;
        end
    end
end

%% NonLinear Coupler Force Vector
Fosc_NLC =  zeros(DoF_conf, 1); % Positive as a force
for i_NLC = 1:length(alpha_ALL)
    i_OSCs = i_NLC + [0 2];
    dx_NLC = X(i_OSCs(2)) - X(i_OSCs(1));
    F_NLC = alpha_ALL(i_NLC) * dx_NLC.*abs(dx_NLC) ^(beta_ALL(i_NLC)-1);
    
    Fosc_NLC(i_OSCs(1)) = Fosc_NLC(i_OSCs(1)) + F_NLC;
    Fosc_NLC(i_OSCs(2)) = Fosc_NLC(i_OSCs(2)) - F_NLC;
end

%% Compute the Differential Increment
dX = [zeros(DoF_conf), eye(DoF_conf); -M\K, -M\D]*X...
    + [zeros(DoF_conf, 1); M\Fapp_t]...
    + [zeros(DoF_conf, 1); M\Fosc_NLC];
end

