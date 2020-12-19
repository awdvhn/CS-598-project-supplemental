function aResModel = aResFUN_Model(Kcoefs, m, x)
k = Kcoefs(1);
alpha = Kcoefs(2);
beta = Kcoefs(3);

aResModel = -(k/m.*x + alpha/m.*x.*abs(x).^(beta - 1));
end