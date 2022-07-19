%% sensitivity_analysis2
% x = [Kx,Kz,S,P,C,alphas,M,delta_Q];
function alphas=getAlphas(delta_Q)
C = 0.01;
ldm = 0.25;
P = 2;
Q = 3000;
alphas = zeros(1,5);
Qp = [Q delta_Q];
for i = 1:5
    if(i>=2)
        Q=Q+Qp(i);
    end
    part1 = (ldm*Q)^(-P)/C;
    part2 = Qp(i)*sign(Qp(i));
    alphas(i) = log(part2)/log(part1);
end
