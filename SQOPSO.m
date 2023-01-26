
function [Best_Pos,Best_fitness,IterCurve] = SQOPSO(pop,maxIter,lb,ub,dim,fobj)
%% ���ò���
c1 = 2.0;
c2 = 2.0;
w = 0.9;
Vmax = 5;
Vmin = -5;
rand('state',sum(100*clock));% �������ظ������
Ub = ub;
Lb = lb;

k = 0.75; %��ʼ�׶�refrOBL��kֵ
n = 4/3; % ������
det = pi/6;% ������ˮ�ĽǶ� 30��
jmax = pi/2;%����� ���Եݼ�
jmin = 0;

a = 0.1;% ��ż���Ե���ֵ
A = zeros(1,maxIter);% ������
x_m = zeros(1,dim);% ÿһά�ľ�ֵ

%% ��ʼ����Ⱥ�ٶ�|λ��|������Ӧ��ֵ
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
Range = ones(pop,1)*(ub-lb);
X = rand(pop,dim).*Range + ones(pop,1)*lb;    % ��ʼ������Ⱥ
V = rand(pop,dim)*(Vmax-Vmin) + Vmin;         % ��ʼ���ٶ�
fun = fobj; %��Ӧ�Ⱥ���
fitness = zeros(1,pop);
for i = 1:pop
    fitness(i) = fun(X(i,:));
end

%% ��ʼ��OBL
FD = 0;%������
dmc = DMC(X,fitness,pop);
for i = 1:pop
    %���㷴���
%     X_r(i,:) = (Ub + Lb)/2 +(Ub + Lb)/(2*k*n) - X(i,:)./(k*n);
    X_r(i,:) = dmc + dmc/(k*n) - X(i,:)/(k*n);
    %���߽�Լ��
    X_r(i,:) = BoundaryCheck(X_r(i,:),ub,lb,dim);
    %������������Ӧ��ֵ
    fitness_r(i) = fun(X_r(i,:));
    %�Ƚ��滻
    if fitness(i) > fitness_r(i)%��Ӧ��ԽСԽ��
        X(i) = X_r(i);
        fitness(i) = fitness_r(i);
        FD = FD + 1;
    end
end
fprintf('��ʼ�׶θ�����%d������\n',FD);

%% ����ʼ��Ⱥ��Ϊ��ʷ����
pBest = X;
pBestFitness = fitness;

%% ��¼��ʼȫ�����Ž�,Ĭ���Ż���Сֵ��
%Ѱ����Ӧ����С��λ��
[~,index] = min(fitness);
%��¼��Ӧ��ֵ��λ��
gBestFitness = fitness(index);
gBest = X(index,:);
Xnew = X; %��λ��
fitnessNew = fitness;%��λ����Ӧ��ֵ
IterCurve = zeros(1,maxIter);
Iter = zeros(1,maxIter);
LD = 0;%������

%% ��ʼ����
for t = 1:maxIter
    % ���������i����������kȡֵ--���Եݼ�
    j = jmax - (jmax - jmin)*t/maxIter;
    k = 1.4 - 1.4*t/maxIter;
    w = 0.9 - 0.5*t/maxIter;
    %% ���ڶ����ԵĽ׶β���ѡ��
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
    %���¶�̬�������dmc
    dmc = DMC(X,fitness,pop);
    for i = 1:pop
        %�ٶȸ���
        r1 = rand(1,dim);
        r2 = rand(1,dim);
        V(i,:) = w.*V(i,:) + c1.*r1.*(pBest(i,:) - X(i,:)) + c2.*r2.*(gBest - X(i,:));
        %�ٶȱ߽��鼰Լ��
        V(i,:) = BoundaryCheck(V(i,:),Vmax,Vmin,dim);
        %λ�ø���
        Xnew(i,:) = X(i,:) + V(i,:);
        if (Iter(t)/Imax) < a % ̽���̶ȴ�����ֵa������г�������
            X_r(i,:) = -(X(i,:)*tan(det)*(n^2-1)*(tan(j)^3)) - X(i,:);
        else % ����˫�ؾ�ֵ���ĵ�׼����
            %���������
            X_r(i,:) = dmc + dmc/(k*n) - X(i,:)/(k*n);
        end
        %������Ӧ��ֵ
        fitness_r(i) = fun(X_r(i,:));
        fitnessNew(i) = fun(Xnew(i,:));
        %�Ƚ��滻
        if fitnessNew(i) > fitness_r(i)
            Xnew(i) = X_r(i);
            fitnessNew(i) = fitness_r(i);
            LD = LD + 1;
        end
        %λ�ñ߽��鼰Լ��
        Xnew(i,:) = BoundaryCheck(Xnew(i,:),ub,lb,dim);
        fitnessNew(i) = fun(Xnew(i,:));
        %������ʷ����ֵ
        if fitnessNew(i) < pBestFitness(i)
            pBest(i,:) = Xnew(i,:);
            pBestFitness(i) = fitnessNew(i);
        end
        %����ȫ������ֵ
        if fitnessNew(i) < gBestFitness
            gBestFitness = fitnessNew(i);
            gBest = Xnew(i,:);
        end
    end
    X = Xnew;
    fitness = fitnessNew;
    %% ��¼��ǰ��������ֵ��������Ӧ��ֵ
    %��¼���Ž�
    Best_Pos = gBest;
    %��¼���Ž����Ӧ��ֵ
    Best_fitness = gBestFitness;
    %��¼��ǰ���������Ž���Ӧ��ֵ
    IterCurve(t) = gBestFitness;
end
fprintf('�Ŷ��׶θ�����%d������\n',LD);
end