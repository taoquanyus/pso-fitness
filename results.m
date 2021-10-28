format long;
load('t.mat');                               %%%%%��ˮʱ�̣�ʱ��s%%%%%
load('sc.mat');                              %%%%%�۲⽵���λm%%%%%
iteration_times=50;                          %����Ҷ��͵Ĵ���
Q=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];    %%%%%��ˮ����m?/s%%%%%
Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4)];
Time=[0,7200,14400,21600,28800];             %%%%%���ݽڵ㣬ʱ��s%%%%%
z=65;
M=74;                                        %%%%%��ˮ����%%%%%
l=72;                                        %%%%%��ˮ���²����%%%%%
d=60;                                        %%%%%��ˮ���ϲ����%%%%%
r=0.105;                                     %%%%%��װ뾶%%%%%
sa=zeros(1,length(t));
sl=zeros(1,length(t));
sz=zeros(1,length(z));
s=zeros(1,length(t));
% x=[0.000585764066514 0.000123132275533 0.005415238111098 1.408037862239064 0.414684292768394 -2.43 -1.58 -1.58 -1.65 -1.71];
% x=[0.000361023499573133 0.000361023499573133 0.00629988342647045 1.87491379721331 0.769141667845900 0.265557117968485 0.944153188804980 1.66172651960770 2.08151576197445 1.77759318925459];

% x=[0.000376063370815
%    0.000376063370815
%    0.006795013194178
%    1.642323938199680
%    0.805822283487601
%    0.053943766909617
%    0.035954813115569
%   -1.057225576225253
%    1.539969165216566
%   -1.425306512942503]
x=[0.000369227810026   0.000635709590273   0.002740057468125   1.025868560264854   0.233309634775674  -1.672249673077384   2.370211494655230   2.720442825170914 2.613709272743112   1.096723315438424];
P=x(4);
C=x(5);
alpha=[x(6),x(7),x(8),x(9),x(10)];
for p=1:length(t)
    T=x(1)*M;
    a=T/x(3);
    
    if(p>0&&p<20)                               %%%%�����������%%%%%
        Tp=t(p)-Time(1);                         %%%%�ж�ʱ��ڵ�%%%%%
        ur=r^2/(4*a*Tp);
        fun1 = @(y) exp(-y)./y;
        W_ur= integral(fun1,ur,Inf);
        
        for m=1:length(z)
            Sum1=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%��������̩˹������%%%%%
                W_2=integral(fun2,ur,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%����Ҷ��ͺ���%%%%%
                Sum1=Sum1+Fourier*W_2;
            end
            Sum1=Sum1*2*M/(pi*(l-d));
            sa(p)=(Qp(1)/(4*pi*T))*(W_ur+Sum1);           %%%%��һ���ݽ׶ν���%%%%%
            sz(m)=sa(p);
        end
        
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
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%��������̩˹������%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%����Ҷ��ͺ���%%%%%
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            sa(p)=Sum1+Sum2;     %%%%�ڶ����ݽ׶ν���+��һ���ݽ���%%%%%
            sz(m)=sa(p);
        end

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
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%��������̩˹������%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%����Ҷ��ͺ���%%%%%
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
            sa(p)=Sum1+Sum2+Sum3;     %%%%�ڶ����ݽ׶ν���+��һ���ݽ���%%%%%
            sz(m)=sa(p);
        end
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
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%��������̩˹������%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%����Ҷ��ͺ���%%%%%
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
            sa(p)=Sum1+Sum2+Sum3+Sum4;     %%%%�ڶ����ݽ׶ν���+��һ���ݽ���%%%%%
            sz(m)=sa(p);
        end
     
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
                fun2= @(y) (1./y).*(exp(-y-(x(2)/x(1)).*(n*pi*r/M)^2./(4.*y))); %%%%%��������̩˹������%%%%%
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                W_2_5=integral(fun2,ur5,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); %%%%����Ҷ��ͺ���%%%%%
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
            sa(p)=Sum1+Sum2+Sum3+Sum4+Sum5;     %%%%�ڶ����ݽ׶ν���+��һ���ݽ���%%%%%
            sz(m)=sa(p);
        end
        
        B=(Qp(5)*sign(Qp(5)))^alpha(5)*((l-d)/M)^P*C;
        sl(p)=B*(Q(5)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    end
end
F=(1/length(t)*sum((s-sc').^2))^(1/2);
plot(t,sc,'o',t,s,'-');
