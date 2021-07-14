format long;
load('t.mat');                               %%%%%抽水时刻，时间s%%%%%
load('sc.mat');                              %%%%%观测降深，单位m%%%%%
iteration_times=50;                          %傅立叶求和的次数
Q=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];    %%%%%抽水流量m?/s%%%%%
Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4)];
Ql=[-Q(1) Q(1)-Q(2) Q(2)-Q(3) Q(3)-Q(4) Q(4)-Q(5)];
Time=[0,7200,14400,21600,28800];             %%%%%阶梯节点，时间s%%%%%
z=65;
%z=[60 62 64 66 68 70];                      %%%%z代表垂向位置%%%%%%
M=74;                                        %%%%%含水层厚度%%%%%
l=72;                                        %%%%%滤水管下部深度%%%%%
d=60;                                        %%%%%滤水管上部深度%%%%%
r=0.105;                                     %%%%%钻孔半径%%%%%
sa=zeros(1,length(t));
sl=zeros(1,length(t));
sz=zeros(1,length(z));
s=zeros(1,length(t));
% x=[0.000585764066514 0.000123132275533 0.005415238111098 1.408037862239064 0.414684292768394 -2.43 -1.58 -1.57 -1.65 -1.71];
% x=[0.000377504131960592,0.000242980319924702,0.00703323301883215,2.21419888060706,0.284751853811488,-2.43 -1.58 -1.58 -1.65 -1.71];
x=[0.000585764066514 0.000123132275533 0.005415238111098 0 0 0 0 0 0 0];
P=x(4);
C=x(5);
alpha=[x(6),x(7),x(8),x(9),x(10)];
for p=1:length(t)
    T=x(1)*M;
    a=T/x(3);
    %                 ur=r^2/(4*a*t(p));
    %                 fun1 = @(y) exp(-y)./y;      %%%%%泰斯井函数%%%%%
    %                 W_ur= integral(fun1,ur,Inf);
    %                 sum=0;                       %%%%傅里叶求和函数初始值%%%%%
    %                 for m=1:length(z)            %%%%将井筛离散化%%%%%
    %                 for n=1:100
    %                     fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
    %                     W_2=integral(fun2,ur,Inf);
    %                     Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
    %                     sum=sum+Fourier*W_2;
    %                 end
    %                 sum=sum*2*M/pi*(l-d);
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc;
% tic;
% format long;
% load('t.mat');                               %%%%%抽水时刻，时间s%%%%%
% load('sc.mat');                              %%%%%观测降深，单位m%%%%%
% iteration_times=50;                          %傅立叶求和的次数
% Q=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];    %%%%%抽水流量m?/s%%%%%
% Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4)];
% Time=[0,7200,14400,21600,28800];             %%%%%阶梯节点，时间s%%%%%
% %z=65;
% z=[60 62 64 66 68 70];                       %%%%z代表垂向位置%%%%%%
% M=74;                                        %%%%%含水层厚度%%%%%
% l=72;                                        %%%%%滤水管下部深度%%%%%
% d=60;                                        %%%%%滤水管上部深度%%%%%
% r=0.105;                                     %%%%%钻孔半径%%%%%
% s=zeros(1,length(t));
% sz=zeros(1,length(z));
% %x=[ 0.000469030531588 0.000410078248507 0.003264115676654];
% x=[ 0.00033105 0.00021739 0.002];
% for p=1:length(t)
%     T=x(1)*M;
%     a=T/x(3);
%     %                 ur=r^2/(4*a*t(p));
%     %                 fun1 = @(y) exp(-y)./y;      %%%%%泰斯井函数%%%%%
%     %                 W_ur= integral(fun1,ur,Inf);
%     %                 sum=0;                       %%%%傅里叶求和函数初始值%%%%%
%     %                 for m=1:length(z)            %%%%将井筛离散化%%%%%
%     %                 for n=1:100
%     %                     fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%     %                     W_2=integral(fun2,ur,Inf);
%     %                     Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%     %                     sum=sum+Fourier*W_2;
%     %                 end
%     %                 sum=sum*2*M/pi*(l-d);
%     
%     if(p>0&&p<20)                               %%%%定义阶梯流量%%%%%
%         Tp=t(p)-Time(1);                         %%%%判断时间节点%%%%%
%         ur=r^2/(4*a*Tp);
%         fun1 = @(y) exp(-y)./y;
%         W_ur= integral(fun1,ur,Inf);
%         
%         for m=1:length(z)
%             Sum1=0;
%             for n=1:iteration_times
%                 fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%                 W_2=integral(fun2,ur,Inf);
%                 Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%                 Sum1=Sum1+Fourier*W_2;
%             end
%             Sum1=Sum1*2*M/(pi*(l-d));
%             s(p)=(Qp(1)/(4*pi*T))*(W_ur+Sum1);           %%%%第一阶梯阶段降深%%%%%
%             sz(m)=s(p);
%         end
%         s(p)=mean(sz);
%     elseif(p>=20&&p<39)
%         Tp2=t(p)-Time(2);
%         ur1=r^2/(4*a*t(p));
%         ur2=r^2/(4*a*Tp2);
%         
%         fun1 = @(y) exp(-y)./y;
%         W_ur1= integral(fun1,ur1,Inf);
%         W_ur2= integral(fun1,ur2,Inf);
%         
%         for m=1:length(z)
%             Sum1=0;
%             Sum2=0;
%             for n=1:iteration_times
%                 fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%                 W_2_1=integral(fun2,ur1,Inf);
%                 W_2_2=integral(fun2,ur2,Inf);
%                 Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%                 Sum1=Sum1+Fourier*W_2_1;
%                 Sum2=Sum2+Fourier*W_2_2;
%             end
%             
%             Sum1=Sum1*2*M/(pi*(l-d));
%             Sum2=Sum2*2*M/(pi*(l-d));
%             Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
%             Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
%             s(p)=Sum1+Sum2;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
%             sz(m)=s(p);
%         end
%         s(p)=mean(sz);
%     elseif(p>=39&&p<58)
%         Tp2=t(p)-Time(2);
%         Tp3=t(p)-Time(3);
%         ur1=r^2/(4*a*t(p));
%         ur2=r^2/(4*a*Tp2);
%         ur3=r^2/(4*a*Tp3);
%         fun1 = @(y) exp(-y)./y;
%         W_ur1= integral(fun1,ur1,Inf);
%         W_ur2= integral(fun1,ur2,Inf);
%         W_ur3= integral(fun1,ur3,Inf);
%         
%         for m=1:length(z)
%             Sum1=0;
%             Sum2=0;
%             Sum3=0;
%             for n=1:iteration_times
%                 fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%                 W_2_1=integral(fun2,ur1,Inf);
%                 W_2_2=integral(fun2,ur2,Inf);
%                 W_2_3=integral(fun2,ur3,Inf);
%                 Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%                 Sum1=Sum1+Fourier*W_2_1;
%                 Sum2=Sum2+Fourier*W_2_2;
%                 Sum3=Sum3+Fourier*W_2_3;
%             end
%             
%             Sum1=Sum1*2*M/(pi*(l-d));
%             Sum2=Sum2*2*M/(pi*(l-d));
%             Sum3=Sum3*2*M/(pi*(l-d));
%             Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
%             Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
%             Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
%             s(p)=Sum1+Sum2+Sum3;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
%             sz(m)=s(p);
%         end
%         s(p)=mean(sz);
%     elseif(p>=58&&p<77)
%         Tp2=t(p)-Time(2);
%         Tp3=t(p)-Time(3);
%         Tp4=t(p)-Time(4);
%         ur1=r^2/(4*a*t(p));
%         ur2=r^2/(4*a*Tp2);
%         ur3=r^2/(4*a*Tp3);
%         ur4=r^2/(4*a*Tp4);
%         fun1 = @(y) exp(-y)./y;
%         W_ur1= integral(fun1,ur1,Inf);
%         W_ur2= integral(fun1,ur2,Inf);
%         W_ur3= integral(fun1,ur3,Inf);
%         W_ur4= integral(fun1,ur4,Inf);
%         
%         for m=1:length(z)
%             Sum1=0;
%             Sum2=0;
%             Sum3=0;
%             Sum4=0;
%             for n=1:iteration_times
%                 fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%                 W_2_1=integral(fun2,ur1,Inf);
%                 W_2_2=integral(fun2,ur2,Inf);
%                 W_2_3=integral(fun2,ur3,Inf);
%                 W_2_4=integral(fun2,ur4,Inf);
%                 Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%                 Sum1=Sum1+Fourier*W_2_1;
%                 Sum2=Sum2+Fourier*W_2_2;
%                 Sum3=Sum3+Fourier*W_2_3;
%                 Sum4=Sum4+Fourier*W_2_4;
%             end
%             
%             Sum1=Sum1*2*M/(pi*(l-d));
%             Sum2=Sum2*2*M/(pi*(l-d));
%             Sum3=Sum3*2*M/(pi*(l-d));
%             Sum4=Sum4*2*M/(pi*(l-d));
%             Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
%             Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
%             Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
%             Sum4=(Qp(4)/(4*pi*T))*(W_ur4+Sum4);
%             s(p)=Sum1+Sum2+Sum3+Sum4;     
%             sz(m)=s(p);
%         end
%         s(p)=mean(sz);
%     elseif(p>=77)
%         Tp2=t(p)-Time(2);
%         Tp3=t(p)-Time(3);
%         Tp4=t(p)-Time(4);
%         Tp5=t(p)-Time(5);
%         ur1=r^2/(4*a*t(p));
%         ur2=r^2/(4*a*Tp2);
%         ur3=r^2/(4*a*Tp3);
%         ur4=r^2/(4*a*Tp4);
%         ur5=r^2/(4*a*Tp5);
%         fun1 = @(y) exp(-y)./y;
%         W_ur1= integral(fun1,ur1,Inf);
%         W_ur2= integral(fun1,ur2,Inf);
%         W_ur3= integral(fun1,ur3,Inf);
%         W_ur4= integral(fun1,ur4,Inf);
%         W_ur5= integral(fun1,ur5,Inf);
%         
%         for m=1:length(z)
%             Sum1=0;
%             Sum2=0;
%             Sum3=0;
%             Sum4=0;
%             Sum5=0;
%             for n=1:iteration_times
%                 fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%非完整井泰斯井函数%%%%%
%                 W_2_1=integral(fun2,ur1,Inf);
%                 W_2_2=integral(fun2,ur2,Inf);
%                 W_2_3=integral(fun2,ur3,Inf);
%                 W_2_4=integral(fun2,ur4,Inf);
%                 W_2_5=integral(fun2,ur5,Inf);
%                 Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%傅里叶求和函数%%%%%
%                 Sum1=Sum1+Fourier*W_2_1;
%                 Sum2=Sum2+Fourier*W_2_2;
%                 Sum3=Sum3+Fourier*W_2_3;
%                 Sum4=Sum4+Fourier*W_2_4;
%                 Sum5=Sum5+Fourier*W_2_5;
%             end
%             
%             Sum1=Sum1*2*M/(pi*(l-d));
%             Sum2=Sum2*2*M/(pi*(l-d));
%             Sum3=Sum3*2*M/(pi*(l-d));
%             Sum4=Sum4*2*M/(pi*(l-d));
%             Sum5=Sum5*2*M/(pi*(l-d));
%             Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
%             Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
%             Sum3=(Qp(3)/(4*pi*T))*(W_ur3+Sum3);
%             Sum4=(Qp(4)/(4*pi*T))*(W_ur4+Sum4);
%             Sum5=(Qp(5)/(4*pi*T))*(W_ur5+Sum5);
%             s(p)=Sum1+Sum2+Sum3+Sum4+Sum5;     
%             sz(m)=s(p);
%         end
%         s(p)=mean(sz);
%     end
% end
% F=(1/length(t)*sum((s-sc').^2))^(1/2);
% plot(t,sc,'o',t,s,'-');
% toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%考虑井损后%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

