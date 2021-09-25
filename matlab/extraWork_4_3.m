%% 非运算模式选择
%Code with Matlab 2019b
clear;
clc;
plot_only = 1;
%% 例6-4[有限时间定常跟踪系统]
%一阶非线性离散系统状态方程： x(k+1) = x(k) + [x(k).^2 + u(k)].*T x(0) = 3
%约束条件 |u(k)|<=2,|x(k)|<=5,k->[0,1]
%性能指标 J = sum(|x(k)-u(k).^3|.*T)
%设定三级决策顺序 J2 -> J1 -> J0   
%得min(J) 求最优控制
if(~plot_only)
    T = 0.1;
    syms x1 x0 x2 J0 J1 J2 u0 u1 u2
    J2 = abs(x2-u2.^3).*T;         %符号函数不适合大规模矩阵运算，利用subs函数会导致带值效率极低，速度极慢。
                                   %所以采用matlabFunction将符号函数转化为内置函数并进行矩阵计算    
                                   %进一步的，由于CPU计算速度较慢，采用GPU加速
    J2 = matlabFunction(J2)
    x2_ve = gpuArray(-10:0.01:10);
    u2_ve = gpuArray(-2:0.01:2);
    J2_values = gpuArray(zeros(size(x2_ve,2),size(u2_ve,2)));  
    disp('Caculating J2_values Matrix');
    for j = 1:size(x2_ve,2)
        temp = [];
        for i = 1:size(u2_ve,2)
            J2_values(j,i) = J2(u2_ve(i),x2_ve(j));
        end
    end
    J2_values = gather(J2_values);
    x2_ve = gather(x2_ve);
    u2_ve = gather(u2_ve);
    
    J1 = abs(x1-u1.^3).*T
    J1 = matlabFunction(J1)
    x1_ve = gpuArray(-10:0.01:10);
    u1_ve = gpuArray(-2:0.01:2);
    J1_values = gpuArray(zeros(size(x1_ve,2),size(u1_ve,2)));
    disp('Caculating J1_values Matrix');
    for j = 1:size(x1_ve,2)
        temp = [];
        for i = 1:size(u1_ve,2)
            J1_values(j,i) = J1(u1_ve(i),x1_ve(j));
        end
    end
    J1_values = gather(J1_values);
    x1_ve = gather(x1_ve);
    u1_ve = gather(u1_ve);
    
    x0_ve = 3;
    u0_ve = -2:1:2;
    x01_ve = x0_ve + (x0_ve.^2 + u0_ve).*T;   %u(0) x(0)取值下对应的 x(1)
    %找x(1) 对应的 u(1) 使得J1最小
    %算力可行的情况下，可尽量增加x(1)和u(1)的取值细粒度，尽量不用插值
    u01_ve = [];                    %x01 对应x(1) u01对应u(1)
    disp('Caculating suitable u(1)');
    for z = x01_ve
        z = 3.7;
        index= find(x1_ve(1,:) == z);
        target_array = J1_values(index,:);
        Jmin = min(target_array);
        index = find(target_array == Jmin);
        u01_ve = [u01_ve,u1_ve(index)];
    end
    x12_ve = x01_ve + (x01_ve.^2 + u01_ve).*T;            
    
    J2_min = [];
    u12_ve = [];                   %x12 对应x(2) u12对应u(2)
    disp('Caculating suitable u(2)');
    for z = x12_ve
        index= find(roundn(x2_ve,-2) == roundn(z,-2));
        target_array = J2_values(index,:);
        Jmin = min(target_array);
        J2_min = [J2_min,Jmin];
        index = find(target_array == min(target_array));
        u12_ve = [u12_ve,u2_ve(index)];
    end
    J1_min = abs(x01_ve-u01_ve.^3).*T + J2_min;
    J0_min = abs(x0_ve-u0_ve.^3).*T + J1_min;
    x3_ve = x12_ve + (x12_ve.^2 + u12_ve).*T;  %x(3)
else
    clc;
    clear;
    load('../Res_Data/4-3(1).mat')
    plot_only = 1;
end
min_index = find(J0_min == min(J0_min));
best_x0 = x0
best_x1 = x01_ve(min_index)
best_x2 = x12_ve(min_index)

best_u0 = u0_ve(min_index)
best_u1 = u12_ve(min_index)

min_J0 = J0_min(min_index)
min_J1 = J1_min(min_index)
min_J2 = J2_min(min_index)