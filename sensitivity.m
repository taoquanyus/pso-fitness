%% sensitivity_analysis2
% x = [Kx,Kz,S,P,C,alphas,M,delta_Q];
function F=sensitivity(x)
Kx = x(1);
Kz = x(2);
S = x(3);
P = x(4);
C = x(5);
alpha = x(6:10);
M = 1/x(11);
delta_Q = x(12);

format long;
load('t.mat');                               
iteration_times=50;                          
% Q=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];    
% Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4)];
% Ql=[-Q(1) Q(1)-Q(2) Q(2)-Q(3) Q(3)-Q(4) Q(4)-Q(5)];
Q(1) = 1800;
for i=2:5
    Q(i)=Q(i-1)+delta_Q;
end
Qp = [Q(1) delta_Q delta_Q delta_Q delta_Q];
Ql = [-Q(1) -delta_Q -delta_Q -delta_Q -delta_Q];
Time=[0,7200,14400,21600,28800];
Time=Time/(24*3600);
z=65;
l=72;                                        
d=60;                                        
r=0.105;                                     
sa=zeros(1,length(t));
sl=zeros(1,length(t));
sz=zeros(1,length(z));
s=zeros(1,length(t));
% P=x(4);
% C=x(5);
% alpha=[x(6),x(7),x(8),x(9),x(10)];
% alpha=[-4.680053227263616  -3.046726625655400  x(6) -3.200091454643404  -3.294421565104963];
%-4.680053227263616  -3.046726625655400  -1.218310929245575 -3.200091454643404  -3.294421565104963
for p=1:length(t)
    T=Kx*M;
    a=T/S;

    if(p>0&&p<20)                               
        Tp=t(p)-Time(1);                         
        ur=r^2/(4*a*Tp);
        fun1 = @(y) exp(-y)./y;
        W_ur= integral(fun1,ur,Inf);
        
        for m=1:length(z)
            Sum1=0;
            for n=1:iteration_times
                fun2= @(y) (1./y).*(exp(-y-(Kz/Kz).*(n*pi*r/M)^2./(4.*y))); 
                W_2=integral(fun2,ur,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); 
                Sum1=Sum1+Fourier*W_2;
            end
            Sum1=Sum1*2*M/(pi*(l-d));
            sa(p)=(Qp(1)/(4*pi*T))*(W_ur+Sum1);           
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
                fun2= @(y) (1./y).*(exp(-y-(Kz/Kz).*(n*pi*r/M)^2./(4.*y))); 
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); 
                Sum1=Sum1+Fourier*W_2_1;
                Sum2=Sum2+Fourier*W_2_2;
            end
            
            Sum1=Sum1*2*M/(pi*(l-d));
            Sum2=Sum2*2*M/(pi*(l-d));
            Sum1=(Qp(1)/(4*pi*T))*(W_ur1+Sum1);
            Sum2=(Qp(2)/(4*pi*T))*(W_ur2+Sum2);
            sa(p)=Sum1+Sum2;     
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
                fun2= @(y) (1./y).*(exp(-y-(Kz/Kz).*(n*pi*r/M)^2./(4.*y))); 
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); 
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
            sa(p)=Sum1+Sum2+Sum3;     
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
                fun2= @(y) (1./y).*(exp(-y-(Kz/Kz).*(n*pi*r/M)^2./(4.*y))); 
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); 
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
            sa(p)=Sum1+Sum2+Sum3+Sum4;     
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
                fun2= @(y) (1./y).*(exp(-y-(Kz/Kz).*(n*pi*r/M)^2./(4.*y))); 
                W_2_1=integral(fun2,ur1,Inf);
                W_2_2=integral(fun2,ur2,Inf);
                W_2_3=integral(fun2,ur3,Inf);
                W_2_4=integral(fun2,ur4,Inf);
                W_2_5=integral(fun2,ur5,Inf);
                Fourier=(1/n)*(sin(n*pi*l/M)-sin(n*pi*d/M))*cos(n*pi*z(m)/M); 
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
            sa(p)=Sum1+Sum2+Sum3+Sum4+Sum5;
            sz(m)=sa(p);
        end
        
        B=(Qp(5)*sign(Qp(5)))^alpha(5)*((l-d)/M)^P*C;
        sl(p)=B*(Q(5)^P);
        sa(p)=mean(sz);
        s(p)=sa(p)+sl(p);
    end
end

F=s;
