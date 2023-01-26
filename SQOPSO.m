
function [Best_Pos,Best_fitness,IterCurve] = SQOPSO(pop,maxIter,lb,ub,dim,fobj)
%% 设置参数
c1 = 2.0;
c2 = 2.0;
w = 0.9;
Vmax = 5;
Vmin = -5;
rand('state',sum(100*clock));% 产生非重复随机数
Ub = ub;
Lb = lb;

k = 0.75; %初始阶段refrOBL的k值
n = 4/3; % 折射率
det = pi/6;% 筷子入水的角度 30°
jmax = pi/2;%入射角 线性递减
jmin = 0;

a = 0.1;% 对偶策略的阈值
A = zeros(1,maxIter);% 多样性
x_m = zeros(1,dim);% 每一维的均值

%% 初始化种群速度|位置|计算适应度值
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
Range = ones(pop,1)*(ub-lb);
X = rand(pop,dim).*Range + ones(pop,1)*lb;    % 初始化粒子群
V = rand(pop,dim)*(Vmax-Vmin) + Vmin;         % 初始化速度
fun = fobj; %适应度函数
fitness = zeros(1,pop);
for i = 1:pop
    fitness(i) = fun(X(i,:));
end

%% 初始化OBL
FD = 0;%计数器
dmc = DMC(X,fitness,pop);
for i = 1:pop
    %计算反向解
%     X_r(i,:) = (Ub + Lb)/2 +(Ub + Lb)/(2*k*n) - X(i,:)./(k*n);
    X_r(i,:) = dmc + dmc/(k*n) - X(i,:)/(k*n);
    %检查边界约束
    X_r(i,:) = BoundaryCheck(X_r(i,:),ub,lb,dim);
    %计算折射解的适应度值
    fitness_r(i) = fun(X_r(i,:));
    %比较替换
    if fitness(i) > fitness_r(i)%适应度越小越优
        X(i) = X_r(i);
        fitness(i) = fitness_r(i);
        FD = FD + 1;
    end
end
fprintf('初始阶段更新了%d个样本\n',FD);

%% 将初始种群作为历史最优
pBest = X;
pBestFitness = fitness;

%% 记录初始全局最优解,默认优化最小值。
%寻找适应度最小的位置
[~,index] = min(fitness);
%记录适应度值和位置
gBestFitness = fitness(index);
gBest = X(index,:);
Xnew = X; %新位置
fitnessNew = fitness;%新位置适应度值
IterCurve = zeros(1,maxIter);
Iter = zeros(1,maxIter);
LD = 0;%计数器

%% 开始迭代
for t = 1:maxIter
    % 更新入射角i和收缩因子k取值--线性递减
    j = jmax - (jmax - jmin)*t/maxIter;
    k = 1.4 - 1.4*t/maxIter;
    w = 0.9 - 0.5*t/maxIter;
    %% 基于多样性的阶段策略选择
    div = 0;
    divp = 0;
    for p = 1:dim
        div_mean = mean(X(:,p));
        for q = 1:pop
            divp = divp + sqrt((X(q,p) - div_mean)^2);
        end
        divp = divp/pop;
        div = div + divp;
    end
    Iter(t) = div/dim;
    Imax = Iter(1);
    %更新动态混合中心dmc
    dmc = DMC(X,fitness,pop);
    for i = 1:pop
        %速度更新
        r1 = rand(1,dim);
        r2 = rand(1,dim);
        V(i,:) = w.*V(i,:) + c1.*r1.*(pBest(i,:) - X(i,:)) + c2.*r2.*(gBest - X(i,:));
        %速度边界检查及约束
        V(i,:) = BoundaryCheck(V(i,:),Vmax,Vmin,dim);
        %位置更新
        Xnew(i,:) = X(i,:) + V(i,:);
        if (Iter(t)/Imax) < a % 探索程度大于阈值a，则进行超对立。
            X_r(i,:) = -(X(i,:)*tan(det)*(n^2-1)*(tan(j)^3)) - X(i,:);
        else % 基于双重均值中心的准对立
            %计算折射解
            X_r(i,:) = dmc + dmc/(k*n) - X(i,:)/(k*n);
        end
        %计算适应度值
        fitness_r(i) = fun(X_r(i,:));
        fitnessNew(i) = fun(Xnew(i,:));
        %比较替换
        if fitnessNew(i) > fitness_r(i)
            Xnew(i) = X_r(i);
            fitnessNew(i) = fitness_r(i);
            LD = LD + 1;
        end
        %位置边界检查及约束
        Xnew(i,:) = BoundaryCheck(Xnew(i,:),ub,lb,dim);
        fitnessNew(i) = fun(Xnew(i,:));
        %更新历史最优值
        if fitnessNew(i) < pBestFitness(i)
            pBest(i,:) = Xnew(i,:);
            pBestFitness(i) = fitnessNew(i);
        end
        %更新全局最优值
        if fitnessNew(i) < gBestFitness
            gBestFitness = fitnessNew(i);
            gBest = Xnew(i,:);
        end
    end
    X = Xnew;
    fitness = fitnessNew;
    %% 记录当前迭代最优值和最优适应度值
    %记录最优解
    Best_Pos = gBest;
    %记录最优解的适应度值
    Best_fitness = gBestFitness;
    %记录当前迭代的最优解适应度值
    IterCurve(t) = gBestFitness;
end
fprintf('扰动阶段更新了%d个样本\n',LD);
end