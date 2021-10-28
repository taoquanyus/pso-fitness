function[xm,fv]=PSOClark(clark,N,c1,c2,w_ini,w_end,M,D)
%% %位置限制
tic;
xMax=[1e-2,1e-2,3,1,3,3,3,3,3,3];
xMin=[1e-9,1e-9,1,1e-14,-3,-3,-3,-3,-3,-3];

%% %速度限制
% v_index=0.15;
% vMax=v_index*(xMax-xMin);
% vMax=0.5*ones(1,10);
% vMin=-vMax;
vMax=0.5;
vMin=-0.5;
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
    results(row)=clark(x_cur(row,:));
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
for times=1:M
    w=(w_ini-w_end)*(M-times)/M+w_end;
    
    for row=1:N
        v(row,:)=w*v(row,:)+c1*rand*(pBest(row,:)-x_cur(row,:))+c2*rand*(gBest-x_cur(row,:));
        %速度限制
        for temp=1:D
            if v(row,temp)>vMax
                v(row,temp)=vMax;
            end
            if v(row,temp)<vMin
                v(row,temp)=vMin;
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
        results(row)=clark(x_cur(row,:));
        if results(row)<pBest_result
            pBest_result=results(row);
            pBest(row,:)=x_cur(row,:);
            if pBest_result<gBest_result
                gBest_result=pBest_result;
                nums=nums+1;
                gBest=pBest(row,:);
                result_collection(nums,:)=gBest;
%                 figure(2)
%                 validate(gBest);
            end
        end
    end
    disp("current: "+times+" total: "+M +" time= "+toc+" opt= "+gBest_result+" Update Times= "+nums);
    %     disp(gBest);
    result(times)=gBest_result;
%     figure(1)
%     plot(result);
%     xlabel("lteration Times");
%     ylabel("MSE");
%     pause(0.5);
    
%     if gBest_result<0.1
%         break;
%     end
end

%%
%结束
save fv gBest_result -ascii;
xm=gBest';
fv=gBest_result;