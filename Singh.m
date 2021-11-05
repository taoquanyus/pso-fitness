function mse=Singh(x)
%X [x1,x2,x3,x4]
%x1 C coefficient of well loss
%x2 n loss exponent
%x3 T transimissivity of aquifer
%x4 S storage coefficient of the aquifer

load('t.mat');                               
load('sc.mat');                              
S=x(4);
T=x(3);
r=0.105;%you may change r due to unit
n=x(2);
Q_data=[0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];%m/s
Q_data=Q_data*60;%Unit conversion

a1=1/(4*pi*T);
a2=(r^2*S)/(4*T);
a3=x(1);

format long;

%sw
t=t./60;
sw=zeros(1,length(t));
for p=1:length(t)
    if(p>0&&p<20)
        Q=Q_data(1);
    elseif(p>=20&&p<39)
        Q=Q_data(2);
    elseif(p>=39&&p<58)
        Q=Q_data(3);
    elseif(p>=58&&p<77)
        Q=Q_data(4);
    else
        Q=Q_data(5);
    end
    sw(p)=a3*(Q)^n;
end

%Sa
sa=zeros(1,length(t));
F=zeros(1,length(t));
delta=zeros(1,length(t));
for m=1:length(t)
    if(m==1)
        t_delta=t(1);
    else
        t_delta=t(m)-t(m-1);
    end
    ur=(a2/(m*t_delta));
    fun1 = @(y) exp(-y)./y;
    F(m)= a1*integral(fun1,ur,Inf);
    if(m==1)
        delta(1)=F(1);
    else
        delta(m)=F(m)-F(m-1);
    end
end
for m=1:length(t)
    for k=1:m
        if(k>0&&k<20)
            Q=Q_data(1);
        elseif(k>=20&&k<39)
            Q=Q_data(2);
        elseif(k>=39&&k<58)
            Q=Q_data(3);
        elseif(k>=58&&k<77)
            Q=Q_data(4);
        else
            Q=Q_data(5);
        end
        sa(m)=sa(m)+Q*delta(m-k+1);
    end
end
s=sa+sw;
mse=(1/length(t)*sum((s-sc').^2))^(1/2);

if isnan(mse)
    mse=Inf;
end

plot(t,sc,'o',t,s,'-');

