%% INITIAL
format long;
load('t.mat');                               %%%%%抽水时刻，时间s%%%%%
load('sc.mat');                              %%%%%观测降深，单位m%%%%%
iteration_times=50;                          %傅立叶求和的次数
l=72;
d=60;
M=74;
Q=[1000 1200 1400 1600 1800];
Q=Q/(24*3600);
%%
%alpha
x1=[10^-4 10^-5 10^-5 2 10^-3 0 0 0 0 0];%%alpha=0对比项
x2=[10^-4 10^-5 10^-5 2 10^-3 1 2 3 4 5];%%alpha递增
x3=[10^-4 10^-5 10^-5 2 10^-3 5 4 3 2 1];%%alpha递减
x4=[10^-4 10^-5 10^-5 2 10^-3 5 4 3 2 5];
pars=[x1;x2;x3;x4];

Size=size(pars);
figure(1)
for i=1:Size(1)
    result=fitness2(pars(i,:),Q,l,d,M);
    plot(result);
    hold on;
end
title("alpha");
%%
%Q
x=[10^-4 10^-5 10^-5 2 10^-3 1 1 1 1 1];
Q1=[1000 1200 1400 1600 1800];%%Q递增
Q2=[1000 800 600 400 200];%%Q递减
Q3=[1000 1200 1400 1200 1000];%%Q无规律    
Q_pars=[Q1;Q2;Q3];
Size=size(Q_pars);
figure(2);
for i=1:Size(1)
    result=fitness2(x,Q_pars(i,:),l,d,M);
    plot(result);
    hold on;
end
title("Q");
%%

%     
% 第三部分，l-d/M
Q=[1000 1200 1400 1600 1800];
x=[10^-4 10^-5 10^-5 2 10^-3 1 1 1 1 1];
% (l-d)/M=[0.1 0.5 1];%%l-d/M占比
% 这里调整l-d/m,只需要调整M即可，因为l=72;d=60;M=74;即M=[13.2,66,132];
M=[13.2,66,132];
figure(3);
for i=1:length(M)
    result=fitness2(x,Q_pars(i,:),l,d,M(i));
    plot(result);
    hold on;
end
title('l-d/M');
