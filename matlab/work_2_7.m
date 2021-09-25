%% ������ģʽѡ��
%Code with Matlab 2019b
clear;
clc;
plot_only = 0;
%% ��2-7(1)
%ϵͳ���̣� Dx(t) = u(t)   x(t0) = x0
%����ָ�꣺ J = 0.5 * cx(tf)^2 + 0.5*int(u(t)^2,t0,tf) 
% tf ���� x(tf)����  c > 0
% ������ܶٺ��� 0.5*u(t)^2 + lambda(t)*u(t)
% �����֪����Ϊ����ָ�� + tf�̶� + ĩ������
% �ɵ����򷽳��顢���Ʒ��̡���ֵ����
if(plot_only==1)
    load('../Res_Data/2-7(1).mat')
    plot_only = 1;
    figure(1)
    plot(t,x_eval,'k',t,u_eval,t,lambda_eval,t,J_best_eval,'LineWidth',2);axis([0 2.01 -0.5 1.5])
    xlabel('t/s')
    ylabel('Value')
    title('2-7(1)')
    legend('x^*(t)','u^*(t)','lambda^*(t)','J^*(t)')
else
    syms c x0 tf t0 lambda(t) x(t) u
    eqn = [diff(lambda)==0,diff(x)==-lambda,lambda(t0)==c*x(tf),x(t0)==x0];
    sol = dsolve(eqn,t);
    sol.x
    sol.lambda
    J_best = 0.5*c*subs(sol.x,t,tf)^2 + 0.5*int(sol.lambda.^2,t,t0,tf)
    %���������⣬��ֵ��ͼ
    t = 0:0.1:2;
    x0 = 1;
    tf = 2;
    t0 = 0;
    c = 1;
    x_eval = double(subs(sol.x));
    lambda_eval = double(subs(sol.lambda));
    lambda_eval = repmat(lambda_eval,1,max(size(t))); 
    u_eval = -lambda_eval;
    J_best_eval = double(subs(J_best));
    J_best_eval = repmat(J_best_eval,1,max(size(t))); 
    figure(1)
    plot(t,x_eval,'k',t,u_eval,t,lambda_eval,t,J_best_eval,'LineWidth',2);axis([0 2.01 -0.5 1.5])
    xlabel('t/s')
    ylabel('Value')
    title('2-7(1)')
    legend('x^*(t)','u^*(t)','lambda^*(t)','J^*(t)')
end
%% ��2-7(2)
if(plot_only==1)
    load('../Res_Data/2-7(2).mat')
    plot_only = 1;
    figure(2)
    plot(t,x_eval,'k',t,u_eval,t,lambda_eval,t,J_best_eval,'LineWidth',2);axis([0 2.01 -0.5 1.5])
    xlabel('t/s')
    ylabel('Value')
    title('2-7(2)')
    legend('x^*(t)','u^*(t)','lambda^*(t)','J^*(t)')
else
    clc;
    clear;
    syms c x0 tf t0 lambda(t) x(t) u(t) Dx(t) Du(t) Dx(t)
    L = 0.5*u.^2 + lambda*u;
    fi = 0.5*c*x(tf).^2;
    lag1 = functionalDerivative(L,x)-diff(functionalDerivative(L,Dx),t);
    lag2 = functionalDerivative(L,u)-diff(functionalDerivative(L,Du),t);
    end_f = functionalDerivative(fi,x(tf));
    eqn = [diff(lambda) == lag1,diff(x) == -lambda,lambda(tf) == end_f,x(t0)==x0];
    sol = dsolve(eqn)
    sol.x
    sol.lambda
    J_best = 0.5*c*subs(sol.x,t,tf)^2 + 0.5*int(sol.lambda.^2,t,t0,tf)
    %���������⣬��ֵ��ͼ
    t = 0:0.1:3;
    x0 = 1;
    tf = 3;
    t0 = 0;
    c = 2;
    x_eval = double(subs(sol.x));
    lambda_eval = double(subs(sol.lambda));
    lambda_eval = repmat(lambda_eval,1,max(size(t))); 
    u_eval = -lambda_eval;
    J_best_eval = double(subs(J_best));
    J_best_eval = repmat(J_best_eval,1,max(size(t))); 
    figure(2)
    plot(t,x_eval,'k',t,u_eval,t,lambda_eval,t,J_best_eval,'LineWidth',2);axis([0 2.01 -0.5 1.5])
    xlabel('t/s')
    ylabel('Value')
    title('2-7(2)')
    legend('x^*(t)','u^*(t)','lambda^*(t)','J^*(t)')
end

%% ��2-7(3)
%���ƺ���������

if(plot_only==1)
    figure(3)
    load('../Res_Data/2-7(3).mat')
    plot_only = 1;
    plot(t,x,'k',t,u,'LineWidth',2);axis([0 2.01 -1.5 1.5]);
    xlabel('t(sec.)');
    ylabel('Value')
    title('2-7(3)')
    legend('x^*(t)','u^*(t)');
else
    syms  tf t c x(t) t0 x0
    %��ĩ��״̬����  ��ĩ��ʱ�����ɵ�����£���������������ָ��
    %Ӧ������ⲻ��
    eqn = [diff(x) == -c*x(tf),x(t0) == x0]
    sol = dsolve(eqn,t)
    
    s1=dsolve('Dx=-lambda,Dlambda=0,x(0)=1,lambda(2)=1');

    s1.x
    s1.lambda

    t=0:0.1:2;
    x=double(subs(s1.x));
    u=-double(subs(s1.lambda));
%     x = repmat(x,1,max(size(t)));
    u = repmat(u,1,max(size(t)));
    figure(3)
    plot(t,x,'k',t,u,'LineWidth',2);axis([0 2.01 -1.5 1.5]);
    xlabel('t(sec.)');
    ylabel('Value')
    title('2-7(3)')
    legend('x^*(t)','u^*(t)');
end