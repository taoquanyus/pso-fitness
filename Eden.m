% function rmse = Eden(x)
clear all
load('time.mat')                              
load('s_clark.mat')  

Q_data=[1306 1693 2423 3261 4094 5019];%m^3/day
% Q_data=Q_data/(24*60);%Unit conversion

sw=zeros(length(Tclark));
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
    sw(p)=(2.03*(10^-3)+4.81*(10^-4)*log10(Tclark(p)))*Q+1.82*(10^-7)*Q^2;
end
rmse=sqrt(mean((sw'-Sclark).^2));
plot(Tclark,Sclark,'o',Tclark,sw,'-');