function dmc = DMC(x,fit,popsize)
%% DMC 动态混合中心
%   popsize：种群规模
%   x：初始种群
%   fit：初始种群对应的适应度值
%   T：阈值-使点成为核心点的密度限定值/密度均值
%   Nbr：邻域半径
% clc;
% clear;
% x = linspace(1,10,100);
% det = pi/6;
% dmc = 5;
% j = pi/6;
%     for i = 1:popsize
%         X_r(i) = tan(det)*abs(dmc - x(i))*(((4/3)^2-1)*(tan(j)^3)) + 2*dmc - x(i);       
%     end
% y = linspace(1,100,100);
% figure(1);
% scatter(x,y);
% hold on
% scatter(X_r,y);
% fit= zeros(1,100);
% popsize = 100;
% % T = 5;
% Q = 0.5;
% Nbr = 0.2;

%% 密度阈值T和邻域半径Nbr自适应

u = var(x,0);%样本方差
u = mean(u);
Nbr = u*(4/(3*popsize))^0.2; % 高斯密度核函数
% fprintf('自适应邻域半径Nbr：%f\n',Nbr);

%% 均值中心
x_sum = 0;
for i = 1:popsize
    x_sum = x(i) + x_sum;
end
MC = x_sum / popsize;
% fprintf('均值中心MC = %f\n',MC);

%% 密度中心
k = zeros(1,popsize);%邻域内数据点个数
p = zeros(1,popsize);%密度指标
ind = 1;
n_sum = 0;
for i = 1:popsize
    for j = 1:popsize
       if abs(x(i) - x(j)) < Nbr  %在邻域范围内
           k(i) = k(i) + 1; 
           dist_D =+ sqrt((x(i) - x(j))^2 + (fit(i) - fit(j))^2);
       end
    end
    D(i) =  dist_D / k(i); %平均距离
end
Dmax = max(D);
Dmin = min(D);
for i = 1:popsize
    if k(i) == 1  %去心邻域的密度
        p(i) = 0;
    else
        p(i) = k(i) + (Dmax - D(i)) / (Dmax - Dmin);
    end
end
T = sum(p)/popsize;
for i = 1:popsize
    if p(i) > T %大于阈值，则将该点放入密度数组中。
        x_n(ind) = x(i);
        ind = ind + 1;
    end
end
for i = 1:ind-1
   n_sum = x_n(i) + n_sum;
end
DC = n_sum / ind;%密度中心
% fprintf('密度中心DC = %f\n',DC);
dmc = (MC + DC)/2;
% fprintf('动态混合中心DMC = %f\n',dmc);
end

