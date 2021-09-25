%% 非运算模式选择
%Code with Matlab 2019b
clear;
clc;
plot_only = 0;
%% 例3-2(1)
%系统方程： Dx1(t) = -x1(t)+u(t)   x1(0) = 1
%          Dx2(t) = x1(t)         x2(0) = 0
%控制受约束 |u(t)| <= 1
% J=x2(1)
% x(tf) 自由
% 构造哈密顿函数 H(x,u,lambda) = lambda1(t)*(-x1(t)+u(t))+lambda2(t)*x1(t)
%由性能指标可知 tf = 1
if(plot_only == 1)
    load('../Res_Data/3-2(1).mat')
    figure(1)
    plot(t,lambda1_eval,t,lambda2_eval,'LineWidth',2);axis([0 1 0 1.5]);
    title('lambda curve at [0,1]');xlabel('t/s');ylabel('lambda value');
    legend('lambda1(t)','lambda2(t)');
    figure(2)
    plot(t,x1_eval,t,x2_eval,'LineWidth',2);axis([0 1 -0.5 1.5]);
    title('x1,x2 curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('x1(t)','x2(t)');
    figure(3)
    plot(t,u_eval,'LineWidth',2);axis([0 1.2 -2 1]);
    title('u curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('u');
    figure(4)
    plot(t,J_eval,t,H_eval,'LineWidth',2);axis([0 1 -1 1]);
    title('J H curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('J','H');
else
    syms t lambda1(t) lambda2(t) x1(t) x2(t) u(t) x1_1 x2_1 
    fi = x2_1;
    %协态方程组 + 边界条件   解微分方程
    eqn = [diff(lambda1) == lambda1 - lambda2,diff(lambda2) == 0, ...
                lambda1(1)==diff(fi,x1_1),lambda2(1)==diff(fi,x2_1)];
    sol = dsolve(eqn);
    sol.lambda1
    sol.lambda2

    %根据极小值原理
    u = -sign(sol.lambda1)

    %得最优控制u*带入系统方程，求得轨线x的解析表达
    eqn1 = [ diff(x1) == -x1 + u,diff(x2) == x1, x1(0) == 1,x2(0) == 0];
    sol1 = dsolve(eqn1)
    sol1.x1
    sol1.x2

    %根据轨线方程x2得到最优性能指标J*=x2(1)    
    J = subs(sol1.x2,1)             
    %根据x1,x2,u,lambda1,lambda2得到最优H*
    H = sol.lambda1*(-sol1.x1+u)+sol.lambda2*sol1.x1

    %带值 画图
    t = 0:0.001:1;
    u_eval = double(subs(u));
    lambda1_eval = double(subs(sol.lambda1));
    lambda2_eval = double(subs(sol.lambda2));
    lambda2_eval = repmat(lambda2_eval,1,max(size(t)));
    x1_eval = double(subs(sol1.x1));
    x2_eval = double(subs(sol1.x2));
    J_eval = double(J);
    J_eval = repmat(J_eval,1,max(size(t)));
    H_eval = double(subs(H));
    figure(1)
    plot(t,lambda1_eval,t,lambda2_eval,'LineWidth',2);axis([0 1 0 1.5]);
    title('lambda curve at [0,1]');xlabel('t/s');ylabel('lambda value');
    legend('lambda1(t)','lambda2(t)');
    figure(2)
    plot(t,x1_eval,t,x2_eval,'LineWidth',2);axis([0 1 -0.5 1.5]);
    title('x1,x2 curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('x1(t)','x2(t)');
    figure(3)
    plot(t,u_eval,'LineWidth',2);axis([0 1.2 -2 1]);
    title('u curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('u');
    figure(4)
    plot(t,J_eval,t,H_eval,'LineWidth',2);axis([0 1 -1 1]);
    title('J H curve at [0,1]');xlabel('t/s');ylabel('Value');
    legend('J','H');
end
