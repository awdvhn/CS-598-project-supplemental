function dX = DynSys_SingleOsc_CubicNonlinear(t, X, Coefs, t_exp, F_exp, m, t_Fend)
k = Coefs(1);
alpha = Coefs(2);
beta = Coefs(3);
d = Coefs(4);

%% Forcing to be Applied
if t < t_Fend
    Fapp_t = interp1(t_exp, F_exp, t);
else
    Fapp_t = 0;
end

%% Compute Derivatives
dX = [0, 1; -k/m, -d/m]*X + [0; 1/m]*Fapp_t...
    +[0; -(alpha/m).*X(1).*abs(X(1)).^(beta-1)];
end

