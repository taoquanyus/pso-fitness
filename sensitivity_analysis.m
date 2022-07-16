%% sensitivity analysis
delta_p = 0.01;
alphas = 0.1*ones(1,5);
delta_Q = 200;
C = 0.01;
P = 2;
Kx = 5;
Kz = 0.5;
S = 0.0001;
M = 72;
rate = 0.01;

x = [Kx,Kz,S,P,C,alphas,M,delta_Q];
load('t.mat');
%% P
delta_x = zeros(1,length(x));
delta_x(4) = x(4)*rate;
partial_C = (sensitivity(x+delta_x)-sensitivity(x))/delta_x(4);
r = x(4)*partial_C;
plot(t,r,'b');
hold on;

%% delta_Q
delta_x = zeros(1,length(x));
delta_x(12) = x(12)*rate;
partial_C = (sensitivity(x+delta_x)-sensitivity(x))/delta_x(12);
r = x(12)*partial_C;
plot(t,r,'g');

%% M/(l-d)
delta_x = zeros(1,length(x));
delta_x(11) = x(11)*rate;
partial_C = (sensitivity(x+delta_x)-sensitivity(x))/delta_x(11);
r = x(11)*partial_C;
plot(t,r,'r');

%% alpha
delta_x = zeros(1,length(x));
delta_x(6:10) = x(6:10)*rate;
partial_C = (sensitivity(x+delta_x)-sensitivity(x))/delta_x(6);
r = x(6)*partial_C;
plot(t,r,'-');