function[x,f]=PSO_Singh(Singh,N,c1,c2,w_ini,w_end,M,D)
tic;
%x1 C coefficient of well loss
%x2 n loss exponent
%x3 T transimissivity of aquifer
%x4 S storage coefficient of the aquifer
xMax=[1e-1,3.5,1000,1e-2];
xMin=[1e-9,1.5,0,1e-9];

vMax=0.5*xMax;
vMin=-0.5*xMin;
x_cur=rand(N,D).*(xMax-xMin)+xMin;
v=rand(N,D).*(vMax-vMin)+vMin;
%% %首次迭代
results=zeros(N,1);%结果矩阵
pBest=x_cur;
gBest=x_cur(1,:);
pBest_result=inf;%局部最优解初始化
gBest_result=inf;%全局最优解初始化

for row=1:N
    results(row)=Singh(x_cur(row,:));
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
nums=1;
result_collection(nums,:)=gBest;

gBest_result=Singh(gBest);
result=zeros(1,M);
for times=1:M
    w=(w_ini-w_end)*(M-times)/M+w_end;
    
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
            x_cur(row,2)=0.1*x_cur(row,1);
        end
     
        flag=0;
        for temp=1:D
            if x_cur(row,temp)>xMax(temp)
                flag=1;
            end
            if x_cur(row,temp)<xMin(temp)
                %                 x_cur(row,temp)=xMin(temp);
                flag=1;
            end
        end
        if(flag==1)
            continue;
        end
        
        %迭代
        results(row)=Singh(x_cur(row,:));
        if results(row)<pBest_result
            pBest_result=results(row);
            pBest(row,:)=x_cur(row,:);
            if pBest_result<gBest_result
                gBest_result=pBest_result;
                nums=nums+1;
                gBest=pBest(row,:);
                result_collection(nums,:)=gBest;
            end
        end
    end
    disp("current: "+times+" total: "+M +" time= "+toc+" opt= "+gBest_result+" Update Times= "+nums);
    %     disp(gBest);
    result(times)=gBest_result;
end
%结束
save fv gBest_result -ascii;
x=gBest';
f=gBest_result;