%% 非运算模式选择
%Code with Matlab 2019b
clear;
clc;
plot_only = 1;
%% 问题5-2(1)[有限时间状态调节器]
%系统方程： Dx(t) = ax(t)+u(t)  x(0) = x0
%性能指标 J = 0.5*f*x(tf)^2+0.5*int([qx(t)^2+ru(t)^2],0,tf);f>=0,q>0,r>0
%得min(J) 求最优控制u(t)
if(~plot_only)
    syms t p(t) a r q f tf u(t) x(t) x0
    eqn = [diff(p)== -2*a*p+(1/r)*p.^2-q,p(tf) == f];          %P166 黎卡迪方程推导
    sol_p = dsolve(eqn,t)
    eqn1 = [diff(x)==(a-(1/r).*sol_p).*x,x(0) == x0]; 
    sol_x = dsolve(eqn1,t)
    u = -(1/r).*sol_p.*sol_x
    %带值画图
    x0 = 1;
    tf = 1;
    a = -1;
    f = 0;
    q = 1;

    t = 0:0.01:tf;
    x_eval = subs(sol_x);
    u_eval = subs(u);
    p_eval =  subs(sol_p);
else
    load('../Res_Data/5-2(1).mat')
end
%% (1.2)绘制r = [1,0.2,0.02,100] 的轨线方程x(t)

figure(1)
for r = [1,0.2,0.02,100]
    temp_x_eval = double(subs(x_eval));
    plot(t,real(temp_x_eval));hold on;
end
fsc = ["r = 100","r = 1","r = 0.2","r = 0.02"];
pos = [[0.6,0.6];[0.52,0.52];[0.42,0.42];[0.30,0.18]];
for j=1:size(fsc,2)
    str = sprintf(fsc(j));
    text(pos(j,1),pos(j,2),str);
end
xlabel('t');ylabel('x(t)');
%% (1.3)绘制r = [1,0.2,0.02,100] 的p(t)
figure(2)
for r = [1,0.2,0.02,100]
    temp_p_eval = double(subs(p_eval));
    plot(t,real(temp_p_eval));hold on;
end
fsc = ["r = 1","r = 0.2","r = 0.02","r = 100"];
pos = [[0.28,0.32];[0.3,0.25];[0.3,0.1];[0.30,0.40]];
for j=1:size(fsc,2)
    str = sprintf(fsc(j));
    text(pos(j,1),pos(j,2),str);
end
xlabel('t');ylabel('p(t)');
%% (1.4)绘制r = [1,0.2,0.02,100] 的u(t)
figure(3)
%绘制r = [1,0.2,0.02,100] 的p(t)
for r = [1,0.2,0.02,100]
    temp_u_eval = double(subs(u_eval));
    plot(t,real(temp_u_eval));hold on;axis([0,1,-5,0.5])
end
fsc = ["r = 1","r = 0.2","r = 0.02","r = 100"];
pos = [[0.1,-0.5];[0.5,-0.6];[0.2,-2.5];[0.20,0.15]];
for j=1:size(fsc,2)
    str = sprintf(fsc(j));
    text(pos(j,1),pos(j,2),str);
end
xlabel('t');ylabel('u(t)');


%% tf = [1,3,5,10] f =[0,1] a = -1 q=r=1
clc;
clearvars -except plot_only
if(~plot_only)
    syms t p(t) a r q f tf u(t) x(t) x0
    eqn = [diff(p)== -2*a*p+(1/r)*p.^2-q,p(tf) == f];          %P166 黎卡迪方程推导
    sol_p = dsolve(eqn,t)
    eqn1 = [diff(x)==(a-(1/r).*sol_p).*x,x(0) == x0]; 
    sol_x = dsolve(eqn1,t)
    u = -(1/r).*sol_p.*sol_x
    %带值画图
    x0 = 1;
    tf_2 = [1,3,5,10];
    f_2 =[0,1];
    a = -1;
    q = 1; r = 1;
    p_eval_2 =  subs(sol_p,{'a','q','r','x0'},{a,q,r,x0});
else
    load('../Res_Data/5-2(2).mat')
end
figure(4)
for z = tf_2
    t = 0:0.01:z;
    p_eval_21 =  subs(p_eval_2,{'t','tf'},{t,z});
    for m = f_2 
        temp_p_eval_21 = double(subs(p_eval_21,{'f'},{m}));
        plot(t,real(temp_p_eval_21));hold on;
    end
end
xlabel('t');ylabel('p(t)');
