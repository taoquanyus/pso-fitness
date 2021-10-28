format long;
load('t.mat');                               %%%%%抽水时刻，时间s%%%%%
load('sc.mat');                              %%%%%观测降深，单位m%%%%%
iteration_times=50;                          %傅立叶求和的次数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第一部分，探讨alpha的影响
Q=[1000 1200 1400 1600 1800];
Q=Q/(24*3600);%将单位设置为m3/s
x(1)=10^-4;
x(2)=10^-5;
x(3)=10^-5;
x(4)=2;
x(5)=10^-3;
%%%开始控制alpha的结果
x=[10^-4 10^-5 10^-5 2 10^-3 0 0 0 0 0];%%alpha=0对比项
x=[10^-4 10^-5 10^-5 2 10^-3 1 2 3 4 5];%%alpha递增
x=[10^-4 10^-5 10^-5 2 10^-3 5 4 3 2 1];%%alpha递减
x=[10^-4 10^-5 10^-5 2 10^-3 5 4 3 2 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第二部分，探讨Q的影响
Q=[1000 1200 1400 1600 1800];%%Q递增
Q=[1000 800 600 400 200];%%Q递减
Q=[1000 1200 1400 1200 1000];%%Q无规律
Q=Q/(24*3600);%将单位设置为m3/s
x(1)=10^-4;
x(2)=10^-5;
x(3)=10^-5;
x(4)=2;
x(5)=10^-3;
%%%开始控制alpha的结果
x=[10^-4 10^-5 10^-5 2 10^-3 1 1 1 1 1];%%alpha=1看形势
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第三部分，探讨l-d/M的影响
Q=[1000 1200 1400 1600 1800];
Q=Q/(24*3600);%将单位设置为m3/s
x(1)=10^-4;
x(2)=10^-5;
x(3)=10^-5;
x(4)=2;
x(5)=10^-3;
%%%开始控制alpha的结果
x=[10^-4 10^-5 10^-5 2 10^-3 1 1 1 1 1];%%alpha=1

(l-d)/M=[0.1 0.5 1];%%l-d/M占比

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];    %%%%%抽水流量m?/s%%%%%
% Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4)];
% Ql=[-Q(1) Q(1)-Q(2) Q(2)-Q(3) Q(3)-Q(4) Q(4)-Q(5)];
Time=[0,7200,14400,21600,28800];             %%%%%阶梯节点，时间s%%%%%
z=65;
% z=[60 62 64 66 68 70 72];                      %%%%z代表垂向位置%%%%%%
% z=[60 61 62 63 64 65 66 67 68 69 70 71 72]; 
M=74;                                        %%%%%含水层厚度%%%%%
l=72;                                        %%%%%滤水管下部深度%%%%%
d=60;                                        %%%%%滤水管上部深度%%%%%
r=0.105;                                     %%%%%钻孔半径%%%%%
sa=zeros(1,length(t));
sl=zeros(1,length(t));
sz=zeros(1,length(z));
s=zeros(1,length(t));
% x=[0.000585764066514 0.000123132275533 0.005415238111098 1.408037862239064 0.414684292768394 -2.43 -1.58 -1.58 -1.65 -1.71];
%x=[0.000511212204736 0.000511212204736 0.000000001000000 1.000000000000000 0.100000000000000 -2.020191292726334 0.391585221065750 0.168130266560995 1.800473137640037 -1.423482035402611];
% x=[0.000511527378169000;0.000511527378169000;1.00000000000000e-09;1.50000000000000;0.100000000000000;-2.75916979531924;-1.54758202711745;1.65022530072580;1.12150128321838;-1.94740711802740];
x=[0.000573956860243   0.000009280131392   0.009992973699525   2.580635919097652   0.037485269053980  -4.680053227263616  -3.046726625655400  -1.218310929245575 -3.200091454643404  -3.294421565104963];
P=x(4);
C=x(5);
alpha=[x(6),x(7),x(8),x(9),x(10)];
for p=1:length(t)
    T=x(1)*M;
    a=T/x(3);
    if(p>0&&p<20)                               %%%%定义阶梯流量%%%%%
        Tp=t(p)-Time(1);                         %%%%判断时间节点%%%%%
        ur=r^2/(4*a*Tp);
        fun1 = @(y) exp(-y)./y;
        W_ur= integral(fun1,ur,Inf);
        
        for m=1:length(z)
            Sum1=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
                W_2=integral(fun2,ur,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
                Sum1=Sum1+Fourier*W_2;
            end
            Sum1=Sum1*2*M/(pi*(l-d));
            sa(p)=(Qp(1)/(4*pi*T))*(W_ur+Sum1);           %%%%第一阶梯阶段降深%%%%%
            sz(m)=sa(p);
        end
        
        %B=(-1)^heaviside(Ql(1))*C*(alpha(1)*(l-d)/M)^P;
        
        %B=Ql(1)^(alpha(1)*heaviside(Ql(1)))*((l-d)/M)^P*C;
        
        %B=alpha(1)^(Qp(1)*sign(Qp(1)))*((l-d)/M)^P*C;
        %B=alpha(1)^(Q(1)*sign(Q(1)))*((l-d)/M)^P*C
        B=(Qp(1)*sign(Qp(1)))^alpha(1)*((l-d)/M)^P*C;
        
        sl(p)=B*(Q(1)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    elseif(p>=20&&p<39)
        Tp2=t(p)-Time(2);
        ur1=r^2/(4*a*t(p));
        ur2=r^2/(4*a*Tp2);
        
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        
        for m=1:length(z)
            Sum1=0;
            Sum2=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            sa(p)=Sum1+Sum2;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
            sz(m)=sa(p);
        end
        %B=(-1)^heaviside(Ql(2))*C*(alpha(2)*(l-d)/M)^P;
        
        %B=alpha(2)^(Qp(2)*sign(Qp(2)))*((l-d)/M)^P*C;
        
        %B=alpha(2)^(Qp(2)*heaviside(Ql(2)))*((l-d)/M)^P*C;
        B=(Qp(2)*sign(Qp(2)))^alpha(2)*((l-d)/M)^P*C;
        sl(p)=B*(Q(2)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    elseif(p>=39&&p<58)
        Tp2=t(p)-Time(2);
        Tp3=t(p)-Time(3);
        ur1=r^2/(4*a*t(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);
        
        for m=1:length(z)
            Sum1=0;
            Sum2=0;
            Sum3=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
                Sum3=Sum3+Fourier*W_2_3;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum3=Sum3*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
            sa(p)=Sum1+Sum2+Sum3;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
            sz(m)=sa(p);
        end
        
        %B=(-1)^heaviside(Ql(3))*C*(alpha(3)*(l-d)/M)^P;
        
        %B=Ql(3)^(alpha(3)*heaviside(Ql(3)))*((l-d)/M)^P*C;
        
        %B=alpha(3)^(Qp(3)*sign(Qp(3)))*((l-d)/M)^P*C;
        
        B=(Qp(3)*sign(Qp(3)))^alpha(3)*((l-d)/M)^P*C;
        sl(p)=B*(Q(3)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    elseif(p>=58&&p<77)
        Tp2=t(p)-Time(2);
        Tp3=t(p)-Time(3);
        Tp4=t(p)-Time(4);
        ur1=r^2/(4*a*t(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        ur4=r^2/(4*a*Tp4);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);
        W_ur4= integral(fun1,ur4,Inf);
        
        for m=1:length(z)
            Sum1=0;
            Sum2=0;
            Sum3=0;
            Sum4=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
                Sum3=Sum3+Fourier*W_2_3;
                Sum4=Sum4+Fourier*W_2_4;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum3=Sum3*2*M/(pi*(l-d));
            Sum4=Sum4*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
            Sum4=(Qp(4)/(4*pi*T))*(W_ur4+Sum4);
            sa(p)=Sum1+Sum2+Sum3+Sum4;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
            sz(m)=sa(p);
        end
        %B=(-1)^heaviside(Ql(4))*C*(alpha(4)*(l-d)/M)^P;
        
        %B=Ql(4)^(alpha(4)*heaviside(Ql(4)))*((l-d)/M)^P*C;
        
        %B=alpha(4)^(Qp(4)*sign(Qp(4)))*((l-d)/M)^P*C;
        
        B=(Qp(4)*sign(Qp(4)))^alpha(4)*((l-d)/M)^P*C;
        sl(p)=B*(Q(4)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    elseif(p>=77)
        Tp2=t(p)-Time(2);
        Tp3=t(p)-Time(3);
        Tp4=t(p)-Time(4);
        Tp5=t(p)-Time(5);
        ur1=r^2/(4*a*t(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        ur4=r^2/(4*a*Tp4);
        ur5=r^2/(4*a*Tp5);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);
        W_ur4= integral(fun1,ur4,Inf);
        W_ur5= integral(fun1,ur5,Inf);
        
        for m=1:length(z)
            Sum1=0;
            Sum2=0;
            Sum3=0;
            Sum4=0;
            Sum5=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                W_2_5=integral(fun2,ur5,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
                Sum3=Sum3+Fourier*W_2_3;
                Sum4=Sum4+Fourier*W_2_4;
                Sum5=Sum5+Fourier*W_2_5;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum3=Sum3*2*M/(pi*(l-d));
            Sum4=Sum4*2*M/(pi*(l-d));
            Sum5=Sum5*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
            Sum4=(Qp(4)/(4*pi*T))*(W_ur4+Sum4);
            Sum5=(Qp(5)/(4*pi*T))*(W_ur5+Sum5);
            sa(p)=Sum1+Sum2+Sum3+Sum4+Sum5;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
            sz(m)=sa(p);
        end
        %B=(-1)^heaviside(Ql(5))*C*(alpha(5)*(l-d)/M)^P;
        
        %B=Ql(5)^(alpha(5)*heaviside(Ql(5)))*((l-d)/M)^P*C;
        
        %B=alpha(5)^(Qp(5)*sign(Qp(5)))*((l-d)/M)^P*C;
        
        B=(Qp(5)*sign(Qp(5)))^alpha(5)*((l-d)/M)^P*C;
        sl(p)=B*(Q(5)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    end
end
F=(1/length(t)*sum((s-sc').^2))^(1/2);
plot(t,sc,'o',t,s,'-');
%         end
%         s(p)=mean(sz);
%     end
% end
% F=(1/length(t)*sum((s-sc').^2))^(1/2);
% plot(t,sc,'o',t,s,'-');
% toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%考虑井损后%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

