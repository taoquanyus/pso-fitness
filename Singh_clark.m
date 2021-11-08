function mse=Singh_clark(x)
%X [x1,x2,x3,x4]
%x1 C coefficient of well loss
%x2 n loss exponent
%x3 T transimissivity of aquifer
%x4 S storage coefficient of the aquifer

load('time.mat')                              
load('s_clark.mat')                            
S=x(4);
T=x(3);
r=0.105;%you may change r due to unit
n=x(2);
Q_data=[1306 1693 2423 3261 4094 5019];%m^3/day
Q_data=Q_data/(24*60);%Unit conversion

a1=1/(4*pi*T);
a2=(r^2*S)/(4*T);
a3=x(1);

format long;

%sw
sw=zeros(1,length(Tclark));
for p=1:length(Tclark)
    if(p>0&&p<27)
        Q=Q_data(1);
    elseif(p>=27&&p<57)
        Q=Q_data(2);
    elseif(p>=57&&p<87)
        Q=Q_data(3);
    elseif(p>=87&&p<117)
        Q=Q_data(4);
    elseif(p>=117&&p<147)
        Q=Q_data(5);
    else
        Q=Q_data(6);
    end
    sw(p)=a3*(Q)^n;
end

%Sa
sa=zeros(1,length(Tclark));
F=zeros(1,length(Tclark));
delta=zeros(1,length(Tclark));
for m=1:length(Tclark)
%     if(m==1)
%         t_delta=Tclark(1);
%     else
%         t_delta=Tclark(m)-Tclark(m-1);
%     end
    t_delta=(1080-5)/175;
    ur=(a2/(m*t_delta));
    fun1 = @(y) exp(-y)./y;
    F(m)= a1*integral(fun1,ur,Inf);
    if(m==1)
        delta(1)=F(1);
    else
        delta(m)=F(m)-F(m-1);
    end
end
for m=1:length(Tclark)
    for k=1:m
        if(k>0&&k<27)
            Q=Q_data(1);
        elseif(k>=27&&k<57)
            Q=Q_data(2);
        elseif(k>=57&&k<87)
            Q=Q_data(3);
        elseif(k>=87&&k<117)
            Q=Q_data(4);
        elseif(k>=117&&k<147)
            Q=Q_data(5);
        else
            Q=Q_data(6);
        end
        sa(m)=sa(m)+Q*delta(m-k+1);
    end
end
s=sa+sw;
mse=(1/length(Tclark)*sum((s-Sclark').^2))^(1/2);

if isnan(mse)
    mse=Inf;
end

plot(Tclark,Sclark,'o',Tclark,s,'-');

