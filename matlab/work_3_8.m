%% 非运算模式选择
%Code with Matlab 2019b
clear;
clc;
plot_only = 1;
%% 问题3-8(1)[时间-燃料最优控制]
%系统方程： Dx1(t) = x2(t)
%          Dx2(t) = u(t) 
%控制受约束 |u(t)| <= 1
%性能指标 J = int(po + |u(t)|,0,tf),其中po为常数且大于0
%初始状态x1(0) = x1_1  x2(0) = x2_1 ， 要求末端状态x1(tf) = 0 x2(tf) = 0
%tf自由  得min(J) 求最优控制u(t)
%构造哈密顿函数 H(x,u,lambda) = po + |u(t)| +lambda1(t)*x2(t)+lambda2(t)*u(t)
%可知u(t) = -dez(lambda2(t))
if(~plot_only)
    syms t lambda1(t) lambda2(t)
    %根据协态方程
    eqn = [diff(lambda1) == 0, diff(lambda2) == -lambda1];
    sol = dsolve(eqn);
    sol.lambda1
    sol.lambda2
end

%H*=0 故lambda1与lambda2 不能同时为0，故控制序列可为:
% u = {+1},{-1},{0,+1},{0,-1},{+1,0,-1},{-1,0,+1}
%开关曲线上取值需要满足 +1:{x1 = 0.5*x2^2}    -1:{x1 = -0.5*x2^2}
%% {+1}情况
if(~plot_only)
    po = 1.6458
    %一步到位 轨线直接从(x1_1,x2_1)转移到(0,0)，设到达时间为tf
    % 0 <= t <= tf,初态(x1_1,x2_1)
    %将u = +1 带入系统方程可得
    %Dx1(t) = x2(t) ,Dx2(t) = 1
    syms tf x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == 1 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_f = subs(sol1.x1,t,tf)
    x2_f = subs(sol1.x2,t,tf)
    sol_x1f = solve(x1_f,tf);
    sol_x2f = solve(x2_f,tf);
    x1_1 = 0.5;
    x2_1 = -1;
    tf_x1f = double(subs(sol_x1f));
    tf_x2f = double(subs(sol_x2f));
    index = find(tf_x1f == tf_x2f);
    if(~isempty(index))
        tf = tf_x1f(index(1))
    else
        disp('incorrect init point!')
    end
    J = po*tf + tf
    t = 0:0.01:tf;
    formatSpec = "Tf = %.2fs";
    str = sprintf(formatSpec,tf);
    formatSpec = "(%.2f , %.2f)";
    str2 = sprintf(formatSpec,x1_1,x2_1);
    figure(1)
    arrowPlot(double(subs(sol1.x1)),double(subs(sol1.x2)),'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('+1 transfer curve');xlabel('x1');ylabel('x2')
    text(0.25,-0.5,str)
    text(x1_1,x2_1-0.1,str2)
    figure(11)
    plot(t,repmat(1,1,max(size(t))),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U +1 curve');xlabel('t/s');

end
%% {-1}情况
if(~plot_only)
    po = 1.6458
    %一步到位 轨线直接从(x1_1,x2_1)转移到(0,0)，设到达时间为tf
    % 0 <= t <= tf,初态(x1_1,x2_1)
    %将u = -1 带入系统方程可得
    %Dx1(t) = x2(t) ,Dx2(t) = -1
    syms tf x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == -1 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_f = subs(sol1.x1,t,tf)
    x2_f = subs(sol1.x2,t,tf)
    sol_x1f = solve(x1_f,tf);
    sol_x2f = solve(x2_f,tf);
    x1_1 = -2;
    x2_1 = 2;
    tf_x1f = double(subs(sol_x1f));
    tf_x2f = double(subs(sol_x2f));
    index = find(tf_x1f == tf_x2f);
    if(~isempty(index))
        tf = tf_x1f(index(1))
    else
        disp('incorrect init point!')
    end
    J = po*tf + tf
    t = 0:0.01:tf;
    formatSpec = "Tf = %.2fs";
    str = sprintf(formatSpec,tf);
    formatSpec = "(%.2f , %.2f)";
    str2 = sprintf(formatSpec,x1_1,x2_1);
    figure(2)
    arrowPlot(double(subs(sol1.x1)),double(subs(sol1.x2)),'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('-1 transfer curve');xlabel('x1');ylabel('x2')
    text(-2,1,str)
    text(x1_1-0.5,x2_1+0.3,str2)
    figure(22)
    plot(t,repmat(-1,1,max(size(t))),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U -1 curve');xlabel('t/s');
end
%% {0,+1}情况
if(~plot_only)
    po = 1.6458
    %两步到位 轨线从(x1_1,x2_1)转移到开关曲线，在到达(0,0)，设到达分别时间为ta,tf.
    %分别将u = {0，+1} 带入系统方程
    % u=0
    syms tf ta x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == 0 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_a = subs(sol1.x1,t,ta)
    x2_a = subs(sol1.x2,t,ta)

    % u=+1
    eqn2 = [diff(x2) == 1 ,diff(x1) == x2,x1(0) == x1_a, x2(0) == x2_a];
    sol2 = dsolve(eqn2,t);
    x1_af = subs(sol2.x1,t,t-ta);
    x2_af = subs(sol2.x2,t,t-ta);
    x1_f = subs(x1_af,t,tf)
    x2_f = subs(x2_af,t,tf)
    sol2_q = solve(x1_f,x2_f,ta,tf);
    x1_1 = 2;
    x2_1 = -1;
    ta = double(subs(sol2_q.ta))
    tf = double(subs(sol2_q.tf))
    J = po * tf - ta 

    t = 0:0.01:tf;
    x1_curve = double(subs(sol1.x1)).*(t<ta) + double(subs(x1_af)).*(t>=ta & t<=tf);
    x2_curve =  repmat(double(subs(sol1.x2)),1,max(size(t))).*(t<ta) + double(subs(x2_af)).*(t>=ta & t<=tf);
    formatSpec = "Ta = %.2fs";
    str = sprintf(formatSpec,ta);
    formatSpec = "Taf = %.2fs";
    str1 = sprintf(formatSpec,tf-ta);
    formatSpec = "(%.1f , %.1f)";
    str2 = sprintf(formatSpec,x1_1,x2_1);
    formatSpec = "(%.1f , %.1f)";
    x1a =double(subs(x1_a));x2a = double(subs(x2_a));
    str3 = sprintf(formatSpec,x1a,x2a);
    figure(3)
    arrowPlot(x1_curve,x2_curve,'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('0 +1 transfer curve');xlabel('x1');ylabel('x2')
    % text(-2,1,str2)
    text(x1_1-0.4,x2_1-0.3,str2)
    text(x1a-0.5,x2a-0.3,str3)
    text((x1_1-x1a)/2 + x1a-0.5,x2a+0.3,str)
    text(-1,-0.5,str1)
    figure(33)
    plot(t,repmat(0,1,max(size(t))).*(t<ta)+repmat(1,1,max(size(t))).*(t<=tf & t>=ta),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [0 +1] curve');xlabel('t/s');
end

%% {0,-1}情况
if(~plot_only)
    po = 1.6458
    %两步到位 轨线从(x1_1,x2_1)转移到开关曲线，在到达(0,0)，设到达分别时间为ta,tf.
    %分别将u = {0，-1} 带入系统方程
    % u=0
    syms tf ta x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == 0 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_a = subs(sol1.x1,t,ta)
    x2_a = subs(sol1.x2,t,ta)

    % u=-1
    eqn2 = [diff(x2) == -1 ,diff(x1) == x2,x1(0) == x1_a, x2(0) == x2_a];
    sol2 = dsolve(eqn2,t);
    x1_af = subs(sol2.x1,t,t-ta);
    x2_af = subs(sol2.x2,t,t-ta);
    x1_f = subs(x1_af,t,tf)
    x2_f = subs(x2_af,t,tf)
    sol2_q = solve(x1_f,x2_f,ta,tf);
    x1_1 = -2.5;
    x2_1 = 1;
    ta = double(subs(sol2_q.ta))
    tf = double(subs(sol2_q.tf))
    J = po * tf - ta 

    t = 0:0.01:tf;
    x1_curve = double(subs(sol1.x1)).*(t<ta) + double(subs(x1_af)).*(t>=ta & t<=tf);
    x2_curve =  repmat(double(subs(sol1.x2)),1,max(size(t))).*(t<ta) + double(subs(x2_af)).*(t>=ta & t<=tf);
    formatSpec = "Ta = %.2fs";
    str = sprintf(formatSpec,ta);
    formatSpec = "Taf = %.2fs";
    str1 = sprintf(formatSpec,tf-ta);
    formatSpec = "(%.1f , %.1f)";
    str2 = sprintf(formatSpec,x1_1,x2_1);
    formatSpec = "(%.1f , %.1f)";
    x1a =double(subs(x1_a));x2a = double(subs(x2_a));
    str3 = sprintf(formatSpec,x1a,x2a);
    figure(4)
    arrowPlot(x1_curve,x2_curve,'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('0 -1 transfer curve');xlabel('x1');ylabel('x2')

    text(x1_1-0.4,x2_1-0.3,str2)
    text(x1a-0.2,x2a+0.3,str3)
    text((x1_1-x1a)/2 + x1a-0.5,x2a+0.3,str)
    text(0.25,0.5,str1)

    figure(44)
    plot(t,repmat(0,1,max(size(t))).*(t<ta)+repmat(-1,1,max(size(t))).*(t<=tf & t>=ta),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [0 -1] curve');xlabel('t/s');
end
%% {-1,0,+1}与{+1,0,-1}情况
if(~plot_only)
    clc;
    clear;
    po = 1.6458

    % {-1,0,+1}
    % u=-1
    syms tf ta tb x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == -1 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_a = subs(sol1.x1,t,ta)
    x2_a = subs(sol1.x2,t,ta)

    % u=0
    eqn2 = [diff(x2) == 0 ,diff(x1) == x2,x1(0) == x1_a, x2(0) == x2_a];
    sol2 = dsolve(eqn2,t);
    x1_ab = subs(sol2.x1,t,t-ta);
    x2_ab = subs(sol2.x2,t,t-ta);
    x1_b = subs(x1_ab,t,tb)
    x2_b = subs(x2_ab,t,tb)

    % u=+1
    eqn3 = [diff(x2) == 1 ,diff(x1) == x2,x1(0) == x1_b, x2(0) == x2_b];
    sol3 = dsolve(eqn3,t);
    x1_bf = subs(sol3.x1,t,t-tb);
    x2_bf = subs(sol3.x2,t,t-tb);
    x1_f = subs(x1_bf,t,tf)
    x2_f = subs(x2_bf,t,tf)


    sol3_q = solve(x1_f,x2_f,x1_a-(po+4)/(2*po)*x2_a.^2,ta,tb,tf);
    x1_1 = 1;
    x2_1 = 1;
    ta = double(subs(sol3_q.ta));
    tb = double(subs(sol3_q.tb));
    tf = double(subs(sol3_q.tf));
    index = find(tf>0);
    if(~isempty(index))
        ta = ta(index)
        tb = tb(index)
        tf = tf(index)
    else
        disp('incorrect init point!')
    end
    J = po * tf + ta + tf - tb 

    t = 0:0.01:tf; 
    x1_curve = double(subs(sol1.x1)).*(t<ta) + double(subs(x1_ab)).*(t>=ta & t<tb) + double(subs(x1_bf)).*(t>=tb & t<=tf);
    x2_curve = double(subs(sol1.x2)).*(t<ta) + double(subs(x2_ab)).*(t>=ta & t<tb) + double(subs(x2_bf)).*(t>=tb & t<=tf);
    formatSpec = "Ta = %.2fs";
    str = sprintf(formatSpec,ta);
    formatSpec = "Tab = %.2fs";
    str1 = sprintf(formatSpec,tb-ta);
    formatSpec = "Tbf = %.2fs";
    str2 = sprintf(formatSpec,tf-tb);
    formatSpec = "(%.1f , %.1f)";
    str3 = sprintf(formatSpec,x1_1,x2_1);
    formatSpec = "(%.1f , %.1f)";
    x1a =double(subs(x1_a));x2a = double(subs(x2_a));
    str4 = sprintf(formatSpec,x1a,x2a);
    formatSpec = "(%.1f , %.1f)";
    x1b =double(subs(x1_b));x2b = double(subs(x2_b));
    str5 = sprintf(formatSpec,x1b,x2b);
    figure(5)
    arrowPlot(x1_curve,x2_curve,'number',2, 'LineWidth', 2,'scale',1.5);axis([-2,2,-2,2]);grid on;hold on;
    text(x1_1-0.3,x2_1+0.1,str3)
    text(x1a-0.3,x1_1-(x2_1-x2a)/2,str)
    text(x1a-0.2,x2a-0.1,str4)
    text((x1_1-x1a)/2 + x1a-0.6,x2a+0.1,str1)
    text(x1b-0.3,x2b-0.1,str5)
    text(x1b-1,x2b/2,str2)


    % {+1,0,-1}
    % u=+1
    po = 1.6458
    syms tf ta tb x1(t) x2(t) u(t) x1_1 x2_1
    eqn1 = [diff(x2) == 1 ,diff(x1) == x2,x1(0) == x1_1, x2(0) == x2_1];
    sol1 = dsolve(eqn1,t);
    x1_a = subs(sol1.x1,t,ta)
    x2_a = subs(sol1.x2,t,ta)

    % u=0
    eqn2 = [diff(x2) == 0 ,diff(x1) == x2,x1(0) == x1_a, x2(0) == x2_a];
    sol2 = dsolve(eqn2,t);
    x1_ab = subs(sol2.x1,t,t-ta);
    x2_ab = subs(sol2.x2,t,t-ta);
    x1_b = subs(x1_ab,t,tb)
    x2_b = subs(x2_ab,t,tb)

    % u=-1
    eqn3 = [diff(x2) == -1 ,diff(x1) == x2,x1(0) == x1_b, x2(0) == x2_b];
    sol3 = dsolve(eqn3,t);
    x1_bf = subs(sol3.x1,t,t-tb);
    x2_bf = subs(sol3.x2,t,t-tb);
    x1_f = subs(x1_bf,t,tf)
    x2_f = subs(x2_bf,t,tf)


    sol3_q = solve(x1_f,x2_f,x1_a+(po+4)/(2*po)*x2_a.^2,ta,tb,tf);
    x1_1 = -1;
    x2_1 = -1;
    ta = double(subs(sol3_q.ta));
    tb = double(subs(sol3_q.tb));
    tf = double(subs(sol3_q.tf));
    index = find(tf>0);
    if(~isempty(index))
        ta = ta(index)
        tb = tb(index)
        tf = tf(index)
    else
        disp('incorrect init point!')
    end
    J = po * tf + ta + tf - tb 

    t = 0:0.01:tf; 
    x1_curve = double(subs(sol1.x1)).*(t<ta) + double(subs(x1_ab)).*(t>=ta & t<tb) + double(subs(x1_bf)).*(t>=tb & t<=tf);
    x2_curve = double(subs(sol1.x2)).*(t<ta) + double(subs(x2_ab)).*(t>=ta & t<tb) + double(subs(x2_bf)).*(t>=tb & t<=tf);
    formatSpec = "Ta = %.2fs";
    str = sprintf(formatSpec,ta);
    formatSpec = "Tab = %.2fs";
    str1 = sprintf(formatSpec,tb-ta);
    formatSpec = "Tbf = %.2fs";
    str2 = sprintf(formatSpec,tf-tb);
    formatSpec = "(%.1f , %.1f)";
    str3 = sprintf(formatSpec,x1_1,x2_1);
    formatSpec = "(%.1f , %.1f)";
    x1a =double(subs(x1_a));x2a = double(subs(x2_a));
    str4 = sprintf(formatSpec,x1a,x2a);
    formatSpec = "(%.1f , %.1f)";
    x1b =double(subs(x1_b));x2b = double(subs(x2_b));
    str5 = sprintf(formatSpec,x1b,x2b);
    arrowPlot(x1_curve,x2_curve,'number',2,'color','r', 'LineWidth', 2,'scale',1.5);axis([-2,2,-2,2]);grid on;hold on;
    title('[-1 0 +1][+1 0 -1] transfer curve');xlabel('x1');ylabel('x2')
    text(x1_1-0.3,x2_1+0.1,str3)
    text(x1a-0.3,x1_1-(x2_1-x2a)/2,str)
    text(x1a-0.2,x2a+0.1,str4)
    text((x1_1-x1a)/2 + x1a,x2a-0.1,str1)
    text(x1b-0.3,x2b+0.1,str5)
    text(x1b+0.3,x2b/2,str2)
        
    x2_z = 0:0.01:2;
    x_zp = -0.5*x2_z.^2;
    x_zn = -((po+4)/(2*po))*x2_z.^2;
    plot(x_zp,x2_z,x_zn,x2_z,'color','k','LineStyle','--');hold on;
    x2_z = -2:0.01:0;
    x_zp = 0.5*x2_z.^2;
    x_zn = ((po+4)/(2*po))*x2_z.^2;
    plot(x_zp,x2_z,x_zn,x2_z,'color','k','LineStyle','--');hold on;

    figure(55)
    plot(t,repmat(-1,1,max(size(t))).*(t<ta)...
          +repmat(0,1,max(size(t))).*(t<tb & t>=ta)...
          +repmat(1,1,max(size(t))).*(t<=tf & t>=tb)...
          ,'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [-1 0 +1] curve');xlabel('t/s');    
    figure(66)
    plot(t,repmat(1,1,max(size(t))).*(t<ta)...
          +repmat(0,1,max(size(t))).*(t<tb & t>=ta)...
          +repmat(-1,1,max(size(t))).*(t<=tf & t>=tb)...
          ,'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [+1 0 -1] curve');xlabel('t/s');
end



%% only plot
if(plot_only)
    clc;
    clear;
    load('../Res_Data/3-8(1_1).mat')
    %u = {+1} condition
    figure(1)
    arrowPlot(double(subs(sol1.x1)),double(subs(sol1.x2)),'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('+1 transfer curve');xlabel('x1');ylabel('x2')
    text(0.25,-0.5,str)
    text(x1_1,x2_1-0.1,str2)
    figure(11)
    plot(t,repmat(1,1,max(size(t))),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U +1 curve');xlabel('t/s');

    
    clc;
    clear;
    %u = {-1} condition
    load('../Res_Data/3-8(1_2).mat')
    figure(2)
    arrowPlot(double(subs(sol1.x1)),double(subs(sol1.x2)),'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('-1 transfer curve');xlabel('x1');ylabel('x2')
    text(-2,1,str)
    text(x1_1-0.5,x2_1+0.3,str2)
    figure(22)
    plot(t,repmat(-1,1,max(size(t))),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U -1 curve');xlabel('t/s');

    clc;
    clear;
    %u = {0,+1} condition
    load('../Res_Data/3-8(1_3).mat')
    figure(3)
    arrowPlot(x1_curve,x2_curve,'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('0 +1 transfer curve');xlabel('x1');ylabel('x2')
    text(x1_1-0.4,x2_1-0.3,str2)
    text(x1a-0.5,x2a-0.3,str3)
    text((x1_1-x1a)/2 + x1a-0.5,x2a+0.3,str)
    text(-1,-0.5,str1)
    figure(33)
    plot(t,repmat(0,1,max(size(t))).*(t<ta)+repmat(1,1,max(size(t))).*(t<=tf & t>=ta),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [0 +1] curve');xlabel('t/s');
    
    clc;
    clear;
    %u = {0,-1} condition
    load('../Res_Data/3-8(1_4).mat')
    figure(4)
    arrowPlot(x1_curve,x2_curve,'number',2,'color', 'b', 'LineWidth', 1,'scale',5);axis([-3,3,-3,3]);grid on;hold on;
    title('0 -1 transfer curve');xlabel('x1');ylabel('x2')
    text(x1_1-0.4,x2_1-0.3,str2)
    text(x1a-0.2,x2a+0.3,str3)
    text((x1_1-x1a)/2 + x1a-0.5,x2a+0.3,str)
    text(0.25,0.5,str1)
    figure(44)
    plot(t,repmat(0,1,max(size(t))).*(t<ta)+repmat(-1,1,max(size(t))).*(t<=tf & t>=ta),'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [0 -1] curve');xlabel('t/s');
    clc;
    clear;
    %u = {0,-1} condition
    load('../Res_Data/3-8(1_5).mat')
    figure(5)
    arrowPlot(x1_curve,x2_curve,'number',2, 'LineWidth', 2,'scale',1.5);axis([-2,2,-2,2]);grid on;hold on;
    text(x1_1-0.3,x2_1+0.1,str3)
    text(x1a-0.3,x1_1-(x2_1-x2a)/2,str)
    text(x1a-0.2,x2a-0.1,str4)
    text((x1_1-x1a)/2 + x1a-0.6,x2a+0.1,str1)
    text(x1b-0.3,x2b-0.1,str5)
    text(x1b-1,x2b/2,str2)    

    
    clc;
    clear;
    %u = {0,-1} condition
    load('../Res_Data/3-8(1_6).mat')
    arrowPlot(x1_curve,x2_curve,'number',2,'color','r', 'LineWidth', 2,'scale',1.5);axis([-2,2,-2,2]);grid on;hold on;
    title('[-1 0 +1][+1 0 -1] transfer curve');xlabel('x1');ylabel('x2')
    text(x1_1-0.3,x2_1+0.1,str3)
    text(x1a-0.3,x1_1-(x2_1-x2a)/2,str)
    text(x1a-0.2,x2a+0.1,str4)
    text((x1_1-x1a)/2 + x1a,x2a-0.1,str1)
    text(x1b-0.3,x2b+0.1,str5)
    text(x1b+0.3,x2b/2,str2)
    x2_z = 0:0.01:2;
    x_zp = -0.5*x2_z.^2;
    x_zn = -((po+4)/(2*po))*x2_z.^2;
    plot(x_zp,x2_z,x_zn,x2_z,'color','k','LineStyle','--');hold on;
    x2_z = -2:0.01:0;
    x_zp = 0.5*x2_z.^2;
    x_zn = ((po+4)/(2*po))*x2_z.^2;
    plot(x_zp,x2_z,x_zn,x2_z,'color','k','LineStyle','--');hold on;
    figure(55)
    plot(t,repmat(-1,1,max(size(t))).*(t<ta)...
          +repmat(0,1,max(size(t))).*(t<tb & t>=ta)...
          +repmat(1,1,max(size(t))).*(t<=tf & t>=tb)...
          ,'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [+1 0 -1] curve');xlabel('t/s');
    figure(66)
    plot(t,repmat(1,1,max(size(t))).*(t<ta)...
          +repmat(0,1,max(size(t))).*(t<tb & t>=ta)...
          +repmat(-1,1,max(size(t))).*(t<=tf & t>=tb)...
          ,'LineWidth', 2);axis([-0.1,(tf+0.1),-2,2]);grid on;
    title('U [+1 0 -1] curve');xlabel('t/s');
end