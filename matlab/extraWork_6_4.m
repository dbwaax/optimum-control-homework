%% 非运算模式选择
%Code with Matlab 2019b
clear;
clc;
plot_only = 1;
%% 例6-4[有限时间定常跟踪系统]
%系统动态方程： Dx(t) = ax(t)+u(t)  x(0) = x0
%              y(t)  = x(t)
%误差方程： e(t) = ym(t)-y(t) = ym(t) - x(t)
%性能指标 J = 0.5*f*e(tf)^2+0.5*int([qx(t)^2+ru(t)^2],0,tf);f>=0,q>0,r>0
%得min(J) 求最优控制u(t)

%% ym(t) = 1(t)
if(~plot_only)
    clc;
    clearvars -except plot_only
    syms a f q r t x(t) u(t) y(t) ym(t) tf g(t) x0 p(t)
    A = a;
    x0 = 0;
    tf = 1;
    q = 1;
    f = 0;
    ym(t) = 1;
    OBSER=obsv(A,1)   %能观性判断
    %能观，存在唯一最优控制u*(t)
    %求解关于P(t)的黎卡迪方程
    eqn = [diff(p) == -2*a*p+(1/r)*p.^2-q, p(tf) == f];
    sol_p = dsolve(eqn,t)
    %求解关于g(t)的微分方程
    eqn1 = [diff(g) == -(a-(1/r).*sol_p).*g - q.*ym, g(tf) == f.*ym(tf)];
    sol_g = dsolve(eqn1,t)
    %求解关于x(t)的微分方程
    eqn2 = [diff(x) == (a-(1/r).*sol_p).*x + (1/r).*sol_g, x(0) == x0];
    sol_x = dsolve(eqn2,t)
    %求u(t)的解析式
    sol_u = -(1/r).*(sol_p.*sol_x-sol_g)
    % a = -1;
    eval_func = [sol_x,sol_g,sol_u];
    a_f = [-1,0,1];
    r_f = [0.01,0.1,1];
    t = 0:0.05:tf;
    temp_list = [];
    count = 1;
    for j = a_f
        for m = eval_func
            figure(count)
            count = count + 1
            for i  = r_f
                temp = double(subs(m,{'a','r','t'},{j,i,t}));
                temp_list = [temp_list;temp];
%                 plot(t,temp);hold on;
            end
        end
    end
else
    clc;
    clear;
    load('../Res_Data/6-4(2).mat')
    plot_only = 1;
    for j=1:3:27
        figure(j)
        plot(t,temp_list(j,:));hold on;
        plot(t,temp_list(j+1,:));hold on;
        plot(t,temp_list(j+2,:));axis([0,1,0,1.5])
    end
       
            
    
end
%% ym(t) = 0.*(t<2.5) + 1.*(t<=5 & t>=2.5) == sign(y=0.4t)
if(~plot_only)
    clc;
    clearvars -except plot_only
    syms a f q r t x(t) u(t) y(t) ym1(t) ym2(t) tf g(t) x0 p(t)
    A = a;
    OBSER=obsv(A,1)   %能观性判断
    tf = 5;
    q = 1;
    f = 0;
    a = 0;
    x0 = 1;
    %能观，存在唯一最优控制u*(t)
    %求解关于P(t)的黎卡迪方程
    options = odeset('RelTol',1e-5,'AbsTol',1e-5);
    [t,p1ode] = ode45(@odefunc2,[tf 0],0,options);
    t = t';p1ode = p1ode';
    ym = 0.*(t<2.5) + 1.*(t>=2.5 & t<=tf);
    %差分方程求解g(k)
    g = zeros(size(p1ode,1),size(p1ode,2));
    g(1) = f.*ym(1);
    r = 1;
    figure(20)
    for j = 1:size(g,2)-1
        g(j+1) =(-(a-(1/r).*p1ode(j)).*g(j) - q.*ym(j)).*(t(j+1)-t(j)) + g(j);
    end
    t = fliplr(t);g = fliplr(g);p1ode = fliplr(p1ode);
    plot(t,g,'LineWidth', 2);hold on;
    options = odeset('RelTol',1e-5,'AbsTol',1e-5);
    [t1,p1ode1] = ode45(@odefunc3,[tf 0],0,options);
    t1 = t1';p1ode1 = p1ode1';
    ym1 = 0.*(t1<2.5) + 1.*(t1>=2.5 & t1<=tf);
    %差分方程求解g(k)
    g1 = zeros(size(p1ode1,1),size(p1ode1,2));
    g1(1) = f.*ym1(1);
    r = 0.01;
    for j = 1:size(g,2)-1
        g1(j+1) =(-(a-(1/r).*p1ode1(j)).*g1(j) - q.*ym1(j)).*(t1(j+1)-t1(j)) + g1(j);
    end
    t1 = fliplr(t1);g1 = fliplr(g1);p1ode1 = fliplr(p1ode1);
    plot(t1,g1,'LineWidth', 2);axis([0,5,-0.1,1.1]);
    xlabel('t/s');ylabel('g(t)');legend('r=1','r=0.01')

    %差分方程求解x(k)
    figure(21)
    x1 = zeros(size(g,1),size(g,2));
    x2 = zeros(size(g1,1),size(g1,2));
    x1(1) = x0;
    x2(1) = x0;
    r = 1;
    for j = 1:size(x1,2)-1
        x1(j+1) = ((a-(1/r).*p1ode(j)).*x1(j)+(1/r).*g(j)).*(t(j+1)-t(j)) + x1(j)        
    end
    plot(t,x1,'LineWidth', 2);hold on;
    r = 0.01;
    for j = 1:size(x2,2)-1
        x2(j+1) = ((a-(1/r).*p1ode1(j)).*x2(j)+(1/r).*g1(j)).*(t1(j+1)-t1(j)) + x2(j)        
    end
    plot(t1,x2,'LineWidth', 2);hold on;
    xlabel('t/s');ylabel('x(t)');legend('r=1','r=0.01')
else
    clc;
    clear;
    load('../Res_Data/6-4(3).mat')
    plot_only = 1;
    figure(223)
    plot(t,g,'LineWidth', 2);hold on;
    plot(t1,g1,'LineWidth', 2);
    xlabel('t/s');ylabel('g(t)');legend('r=1','r=0.01')
    figure(222)
    plot(t,x1,'LineWidth', 2);hold on;
    plot(t1,x2,'LineWidth', 2);
    xlabel('t/s');ylabel('x(t)');legend('r=1','r=0.01')
end
%% ym(t) = 2*sin4t
if(~plot_only)
    clc;
    clearvars -except plot_only
    syms a f q r t x(t) u(t) y(t) ym(t) tf g(t) x0 p(t)
    A = a;
    x0 = 1;
    tf = 5;
    q = 1;
    f = 0;
    a = 0;
    ym(t) = 2*sin(4*t);
    OBSER=obsv(A,1)   %能观性判断
    %能观，存在唯一最优控制u*(t)
    %求解关于P(t)的黎卡迪方程
    eqn = [diff(p) == -2*a*p+(1/r)*p.^2-q, p(tf) == f];
    sol_p = dsolve(eqn,t)
    %求解关于g(t)的微分方程
    eqn1 = [diff(g) == -(a-(1/r).*sol_p).*g - q.*ym, g(tf) == f.*ym(tf)];
    sol_g = dsolve(eqn1,t)
    %求解关于x(t)的微分方程
    eqn2 = [diff(x) == (a-(1/r).*sol_p).*x + (1/r).*sol_g, x(0) == x0];
    sol_x = dsolve(eqn2,t)
    %求u(t)的解析式
    sol_u = -(1/r).*(sol_p.*sol_x-sol_g)
    r_f = [0.01,0.1,1];

    t = 0:0.01:tf;
    tp_list = []
    figure(23)
    for m=r_f
        temp = double(subs(sol_x,{'r','t','a'},{m,t,a}));
        tp_list = [tp_list;temp];
        plot(t,temp,'--');hold on;
    end
    plot(t,double(subs(ym)))
    xlabel('t/s');ylabel('x(t)');legend('r=0.01','r=0.1','r=0.1','y=2sin(4t)')
else
    clc;
    clear;
    load('../Res_Data/6-4(4).mat')
    plot_only = 1;
    figure(24)
    for m=1:size(tp_list,1)
        plot(t,tp_list(m,:),'--');hold on;
    end
    plot(t,double(subs(ym)))
    xlabel('t/s');ylabel('x(t)');legend('r=0.01','r=0.1','r=0.1','y=2sin(4t)')
end