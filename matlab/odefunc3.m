function dp = odefunc2(t,p)
dp = zeros(1,1);
r = 0.01;
q = 1;
a = 0;
dp(1) = -2*a*p(1)+(1/r)*p(1).^2-q;
end