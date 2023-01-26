function dmc = DMC(x,fit,popsize)
%% DMC ��̬�������
%   popsize����Ⱥ��ģ
%   x����ʼ��Ⱥ
%   fit����ʼ��Ⱥ��Ӧ����Ӧ��ֵ
%   T����ֵ-ʹ���Ϊ���ĵ���ܶ��޶�ֵ/�ܶȾ�ֵ
%   Nbr������뾶
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

%% �ܶ���ֵT������뾶Nbr����Ӧ

u = var(x,0);%��������
u = mean(u);
Nbr = u*(4/(3*popsize))^0.2; % ��˹�ܶȺ˺���
% fprintf('����Ӧ����뾶Nbr��%f\n',Nbr);

%% ��ֵ����
x_sum = 0;
for i = 1:popsize
    x_sum = x(i) + x_sum;
end
MC = x_sum / popsize;
% fprintf('��ֵ����MC = %f\n',MC);

%% �ܶ�����
k = zeros(1,popsize);%���������ݵ����
p = zeros(1,popsize);%�ܶ�ָ��
ind = 1;
n_sum = 0;
for i = 1:popsize
    for j = 1:popsize
       if abs(x(i) - x(j)) < Nbr  %������Χ��
           k(i) = k(i) + 1; 
           dist_D =+ sqrt((x(i) - x(j))^2 + (fit(i) - fit(j))^2);
       end
    end
    D(i) =  dist_D / k(i); %ƽ������
end
Dmax = max(D);
Dmin = min(D);
for i = 1:popsize
    if k(i) == 1  %ȥ��������ܶ�
        p(i) = 0;
    else
        p(i) = k(i) + (Dmax - D(i)) / (Dmax - Dmin);
    end
end
T = sum(p)/popsize;
for i = 1:popsize
    if p(i) > T %������ֵ���򽫸õ�����ܶ������С�
        x_n(ind) = x(i);
        ind = ind + 1;
    end
end
for i = 1:ind-1
   n_sum = x_n(i) + n_sum;
end
DC = n_sum / ind;%�ܶ�����
% fprintf('�ܶ�����DC = %f\n',DC);
dmc = (MC + DC)/2;
% fprintf('��̬�������DMC = %f\n',dmc);
end

