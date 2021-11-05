N=1000;%
c1=2;
c2=2;
w_ini=1;
w_end=0.2;
M=2000;
D=4;%变量个数
% [xm,fv]=new_PSO(@fitness,N,c1,c2,w_ini,w_end,M,D);
[x,f]=PSO_Singh(@Singh,N,c1,c2,w_ini,w_end,M,D);
% [xm,fv]=PSOClark(@clark,N,c1,c2,w_ini,w_end,M,D);
x;
f;