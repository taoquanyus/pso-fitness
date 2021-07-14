function[xm,fv]=new_PSO(fitness,N,c1,c2,w,M,D)
%% %位置限制
tic;
xMax=[1e-2,1e-2,1e-2,3,1,3,3,3,3,3];
% xMax=[1e-2,1e-2,1e-2,3,1,3];
xMin=[1e-9,1e-9,1e-9,1,1e-14,-3,-3,-3,-3,-3];
% xMin=[1e-9,1e-9,1e-9,1,1e-14,-3];

%% %速度限制
v_index=0.1;
vMax=v_index*xMax;
vMin=-vMax;

%% %初始化位置,速度

x_cur=rand(N,D).*(xMax-xMin)+xMin;
% x_cur(1,:)=[0.000585764066514 0.000123132275533 0.005415238111098 1.408037862239064 0.414684292768394 -2.43 -1.58 -1.58 -1.65 -1.71];
v=rand(N,D).*(vMax-vMin)+vMin;

%% %首次迭代
results=zeros(N,1);%结果矩阵
pBest=x_cur;
gBest=x_cur(1,:);
pBest_result=inf;%局部最优解初始化
gBest_result=inf;%全局最优解初始化

for row=1:N
    results(row)=fitness(x_cur(row,:));
    if results(row)<pBest_result
        pBest_result=results(row);
        pBest(row,:)=x_cur(row,:);
        
        if pBest_result<gBest_result
            gBest=pBest(row,:);
            gBest_result=pBest_result;
        end
    end
end

%% %M次运算
result=zeros(1,M);
for times=1:M
    for row=1:N
        v(row,:)=w*v(row,:)+c1*rand*(pBest(row,:)-x_cur(row,:))+c2*rand*(gBest-x_cur(row,:));
        
        %速度限制
        for temp=1:D
            if v(row,temp)>vMax(temp)
                v(row,temp)=vMax(temp);
            end
            if v(row,temp)<vMin(temp)
                v(row,temp)=vMin(temp);
            end
        end
        
        %位置限制
        x_cur(row,:)=x_cur(row,:)+v(row,:);
        
        if x_cur(row,2)>x_cur(row,1)
            x_cur(row,2)=x_cur(row,1);
        end
        
        for temp=1:D
            if x_cur(row,temp)>xMax(temp)
                x_cur(row,temp)=xMax(temp);
            end
            if x_cur(row,temp)<xMin(temp)
                x_cur(row,temp)=xMin(temp);
            end
        end
        
        %迭代
        results(row)=fitness(x_cur(row,:));
        if results(row)<pBest_result
            pBest_result=results(row);
            pBest(row,:)=x_cur(row,:);
            if pBest_result<gBest_result
                gBest_result=pBest_result;
                gBest=pBest(row,:);
            end
        end
    end
    disp("current: "+times+" total: "+M +" time= "+toc+" opt= "+gBest_result);
    
    result(times)=gBest_result;
    plot(times,result,'o');
    hold on;
    pause(0.5);
    
    if gBest_result<0.1
        break;
    end
end

%%
%结束
save fv gBest_result -ascii;
xm=gBest';
fv=gBest_result;