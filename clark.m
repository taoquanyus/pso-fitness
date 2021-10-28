function F=clark(x)
format long;
load('time.mat');
load('s_clark.mat');
l=1;
d=0;
M=1;
Q=[1306 1693 2423 3261 4094 5019];    %%%%%抽水流量m?/d%%%%%
Qp=[Q(1) Q(2)-Q(1) Q(3)-Q(2) Q(4)-Q(3) Q(5)-Q(4) Q(6)-Q(5)];
Time=[0,180,360,540,720,900];
r=0.105;                                     %%%%%钻孔半径%%%%%
sa=zeros(1,length(Tclark));
sl=zeros(1,length(Tclark));
s=zeros(1,length(Tclark));
S=x(2);
P=x(3); 
C=x(4);
alpha=[x(5),x(6),x(7),x(8),x(9),x(10)];
for p=1:length(Tclark)
    T=x(1);
    a=T/S;

    if(p>0&&p<27)                               %%%%定义阶梯流量%%%%%
        Tp=Tclark(p)-Time(1);                         %%%%判断时间节点%%%%%
        ur=r^2/(4*a*Tp);
        fun1 = @(y) exp(-y)./y;
        W_ur= integral(fun1,ur,Inf);
        
        sa(p)=(Qp(1)/(4*pi*T))*(W_ur);           %%%%第一阶梯阶段降深%%%%%        
        B=(Qp(1)*sign(Qp(1)))^alpha(1)*((l-d)/M)^P*C;
        sl(p)=B*(Q(1)^P);
        s(p)=sa(p)+sl(p);
    elseif(p>=27&&p<57)
        Tp2=Tclark(p)-Time(2);
        ur1=r^2/(4*a*Tclark(p));
        ur2=r^2/(4*a*Tp2);
        
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        
        Sum1=(Qp(1)/(4*pi*T))*(W_ur1);
        Sum2=(Qp(2)/(4*pi*T))*(W_ur2);
        sa(p)=Sum1+Sum2;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
        B=(Qp(2)*sign(Qp(2)))^alpha(2)*((l-d)/M)^P*C;
        sl(p)=B*(Q(2)^P);
        s(p)=sa(p)+sl(p);
    elseif(p>=57&&p<87)
        Tp2=Tclark(p)-Time(2);
        Tp3=Tclark(p)-Time(3);
        ur1=r^2/(4*a*Tclark(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);

        Sum1=(Qp(1)/(4*pi*T))*(W_ur1);
        Sum2=(Qp(2)/(4*pi*T))*(W_ur2);
        Sum3=(Qp(3)/(4*pi*T))*(W_ur3);
        sa(p)=Sum1+Sum2+Sum3;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
        
        B=(Qp(3)*sign(Qp(3)))^alpha(3)*((l-d)/M)^P*C;
        sl(p)=B*(Q(3)^P);
        s(p)=sa(p)+sl(p);
    elseif(p>=87&&p<117)
        Tp2=Tclark(p)-Time(2);
        Tp3=Tclark(p)-Time(3);
        Tp4=Tclark(p)-Time(4);
        ur1=r^2/(4*a*Tclark(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        ur4=r^2/(4*a*Tp4);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);
        W_ur4= integral(fun1,ur4,Inf);
        
            
        Sum1=(Qp(1)/(4*pi*T))*(W_ur1);
        Sum2=(Qp(2)/(4*pi*T))*(W_ur2);
        Sum3=(Qp(3)/(4*pi*T))*(W_ur3);
        Sum4=(Qp(4)/(4*pi*T))*(W_ur4);
        sa(p)=Sum1+Sum2+Sum3+Sum4;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
        
        B=(Qp(4)*sign(Qp(4)))^alpha(4)*((l-d)/M)^P*C;
        sl(p)=B*(Q(4)^P);
        s(p)=sa(p)+sl(p);
    elseif(p>=117&&p<147)
        Tp2=Tclark(p)-Time(2);
        Tp3=Tclark(p)-Time(3);
        Tp4=Tclark(p)-Time(4);
        Tp5=Tclark(p)-Time(5);
        ur1=r^2/(4*a*Tclark(p));
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
        
        Sum1=(Qp(1)/(4*pi*T))*(W_ur1);
        Sum2=(Qp(2)/(4*pi*T))*(W_ur2);
        Sum3=(Qp(3)/(4*pi*T))*(W_ur3);
        Sum4=(Qp(4)/(4*pi*T))*(W_ur4);
        Sum5=(Qp(5)/(4*pi*T))*(W_ur5);
        sa(p)=Sum1+Sum2+Sum3+Sum4+Sum5;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
        
        B=(Qp(5)*sign(Qp(5)))^alpha(5)*((l-d)/M)^P*C;
        sl(p)=B*(Q(5)^P);
        s(p)=sa(p)+sl(p);
    elseif(p>=147)
        Tp2=Tclark(p)-Time(2);
        Tp3=Tclark(p)-Time(3);
        Tp4=Tclark(p)-Time(4);
        Tp5=Tclark(p)-Time(5);
        Tp6=Tclark(p)-Time(6);
        ur1=r^2/(4*a*Tclark(p));
        ur2=r^2/(4*a*Tp2);
        ur3=r^2/(4*a*Tp3);
        ur4=r^2/(4*a*Tp4);
        ur5=r^2/(4*a*Tp5);
        ur6=r^2/(4*a*Tp6);
        fun1 = @(y) exp(-y)./y;
        W_ur1= integral(fun1,ur1,Inf);
        W_ur2= integral(fun1,ur2,Inf);
        W_ur3= integral(fun1,ur3,Inf);
        W_ur4= integral(fun1,ur4,Inf);
        W_ur5= integral(fun1,ur5,Inf);
        W_ur6= integral(fun1,ur6,Inf);
        
        Sum1=(Qp(1)/(4*pi*T))*(W_ur1);
        Sum2=(Qp(2)/(4*pi*T))*(W_ur2);
        Sum3=(Qp(3)/(4*pi*T))*(W_ur3);
        Sum4=(Qp(4)/(4*pi*T))*(W_ur4);
        Sum5=(Qp(5)/(4*pi*T))*(W_ur5);
        Sum6=(Qp(6)/(4*pi*T))*(W_ur6);
        sa(p)=Sum1+Sum2+Sum3+Sum4+Sum5+Sum6;     %%%%第二阶梯阶段降深+第一阶梯降深%%%%%
        B=(Qp(6)*sign(Qp(6)))^alpha(6)*((l-d)/M)^P*C;
        sl(p)=B*(Q(6)^P);
        s(p)=sa(p)+sl(p);
    end
end
result=(1/length(Tclark)*sum((s-Sclark').^2))^(1/2);
    if isnan(result)
        F=Inf;
    else
        F=result;
    end
end