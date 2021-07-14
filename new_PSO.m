function[xm,fv]=new_PSO(fitness,N,c1,c2,w,M,D)
%% %位置限制
tic;
xMax=[1e-2,1e-2,1e-2,3,1,3,3,3,3,3];
% xMax=[1e-2,1e-2,1e-2,3,1,3];
xMin=[1e-9,1e-9,1e-9,1,1e-14,-3,-3,-3,-3,-3];
% xMin=[1e-9,1e-9,1e-9,1,1e-14,-3];

%% %速度限制
% v_index=0.1;
vMax=1;
vMin=-vMax;

%% %初始化位置,速度

x_cur=rand(N,D).*(xMax-xMin)+xMin;
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
            if v(row,temp)>vMax%(temp)
                v(row,temp)=vMax;%(temp);
            end
            if v(row,temp)<vMin%(temp)
                v(row,temp)=vMin;%(temp);
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
    xlabel=1:length(result);
    plot(xlabel,result);
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




% function[xm,fv]=PSO(fitness,N,c1,c2,w,M,D)
% tic;
% format long;
% % w_ini=0.9;
% % w_end=0.4;
% 
% xMax1=10^-3;   %%%%%%%%%Kr范围，单位m/s%%%%%%%%%
% xMin1=10^-6;
% x1=rand(N,1)*(xMax1-xMin1)+xMin1;
% 
% xMax2=10^-3;    %%%%%%%%%Kz范围,单位m/s%%%%%%%%%
% xMin2=10^-6;
% x2=rand(N,1)*(xMax2-xMin2)+xMin2;
% %Kz应该小于0.1Kr
% 
% xMax3=10^-2;     %%%%%%%%%S范围%%%%%%%%%
% xMin3=10^-6;
% x3=rand(N,1)*(xMax3-xMin3)+xMin3;
% 
% xMax4=3;        %%%%%%%%%P范围%%%%%%%%%
% xMin4=1;
% x4=rand(N,1)*(xMax4-xMin4)+xMin4;
% 
% 
% xMax5=10^-1;     %%%%%%%%%C范围%%%%%%%%%
% xMin5=10^-8;
% x5=rand(N,1)*(xMax5-xMin5)+xMin5;
% 
% 
% xMax6=-2;     %%%%%%%%%alpha1范围%%%%%%%%%
% xMin6=-2.5;
% x6=rand(N,1)*(xMax6-xMin6)+xMin6;
% 
% xMax7=-1.5;     %%%%%%%%%alpha2范围%%%%%%%%%
% xMin7=-2;
% x7=rand(N,1)*(xMax7-xMin7)+xMin7;
% 
% xMax8=-1.5;     %%%%%%%%%alpha3范围%%%%%%%%%
% xMin8=-2.5;
% x8=rand(N,1)*(xMax8-xMin8)+xMin8;
% 
% xMax9=-1.5;     %%%%%%%%%alpha4范围%%%%%%%%%
% xMin9=-2;
% x9=rand(N,1)*(xMax9-xMin9)+xMin9;
% 
% xMax10=-1.5;     %%%%%%%%%alpha5范围%%%%%%%%%
% xMin10=-2;
% x10=rand(N,1)*(xMax10-xMin10)+xMin10;
% 
% % xMax1=20;     %%%%%%%%%Kr范围，单位m/d%%%%%%%%%
% % xMin1=1;
% % x1=rand(N,1)*(xMax1-xMin1)+xMin1;
% % 
% % xMax2=20;    %%%%%%%%%Kz范围,单位m/d%%%%%%%%%
% % xMin2=1;
% % x2=rand(N,1)*(xMax2-xMin2)+xMin2;
% % 
% % xMax3=1e-2;     %%%%%%%%%S范围%%%%%%%%%
% % xMin3=1e-9;
% % x3=rand(N,1)*(xMax3-xMin3)+xMin3;
% % 
% % xMax4=3;        %%%%%%%%%P范围%%%%%%%%%
% % xMin4=1;
% % x4=rand(N,1)*(xMax4-xMin4)+xMin4;
% % 
% % xMax5=1e-1;     %%%%%%%%%C范围%%%%%%%%%
% % xMin5=1e-12;
% % x5=rand(N,1)*(xMax5-xMin5)+xMin5;
% % 
% % xMax6=100;     %%%%%%%%%alpha1范围%%%%%%%%%
% % xMin6=-100;
% % x6=rand(N,1)*(xMax6-xMin6)+xMin6;
% % 
% % xMax7=100;     %%%%%%%%%alpha2范围%%%%%%%%%
% % xMin7=-100;
% % x7=rand(N,1)*(xMax7-xMin7)+xMin7;
% % 
% % xMax8=100;     %%%%%%%%%alpha3范围%%%%%%%%%
% % xMin8=-100;
% % x8=rand(N,1)*(xMax8-xMin8)+xMin8;
% % 
% % xMax9=100;     %%%%%%%%%alpha4范围%%%%%%%%%
% % xMin9=-100;
% % x9=rand(N,1)*(xMax9-xMin9)+xMin9;
% % 
% % xMax10=100;     %%%%%%%%%alpha5范围%%%%%%%%%
% % xMin10=-100;
% % x10=rand(N,1)*(xMax10-xMin10)+xMin10;
% 
% 
% 
% 
% xMax=[xMax1 xMax2 xMax3 xMax4 xMax5 xMax6 xMax7 xMax8 xMax9 xMax10];
% xMin=[xMin1 xMin2 xMin3 xMin4 xMin5 xMin6 xMin7 xMin8 xMin9 xMin10];
% x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10];
% 
% % v_limit=[xMax1/10 xMax2/10 xMax3/10 xMax4/10 xMax5/10 xMax6/10 xMax7/10 xMax8/10 xMax9/10 xMax10/10]; %%速度限制
%  v_limit=xMax;
%  
% %初始速度限制
% v=rand(N,D);
% v(:,1)=v(:,1)*(xMax(1)-xMin(1))+xMin1;
% v(:,2)=v(:,2)*(xMax(2)-xMin(2))+xMin2;
% v(:,3)=v(:,3)*(xMax(3)-xMin(3))+xMin3;
% v(:,4)=v(:,4)*(xMax(4)-xMin(4))+xMin4;
% v(:,5)=v(:,5)*(xMax(5)-xMin(5))+xMin5;
% v(:,6)=v(:,6)*(xMax(6)-xMin(6))+xMin6;
% v(:,7)=v(:,7)*(xMax(7)-xMin(7))+xMin7;
% v(:,8)=v(:,8)*(xMax(8)-xMin(8))+xMin8;
% v(:,9)=v(:,9)*(xMax(9)-xMin(9))+xMin9;
% v(:,10)=v(:,10)*(xMax(10)-xMin(10))+xMin10;
% 
% %%
% %第一次运算
% 
% p=zeros(N,1);
% y=x;
% 
% pg=x(1,:);
% local_opt=inf;
% for i=1:N
%     value=fitness(x(i,:));
%     if isnan(value)
%         p(i)=Inf;%p是处理后的结果矩阵
%     else
%     p(i)=value;
%     y(i,:)=x(i,:);
%         if value<local_opt    %local_opt 全局最优
%             pg=x(i,:);
%             local_opt=value;
%         end
%     end
% end
% 
% %%
% %M次迭代
% result=zeros(1,M);
% for t=1:M
%     for i=1:N
%         v(i,:)=w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
%         if v(i,1)>v_limit(1)
%             v(i,1)=v_limit(1);%%%添加rand*或常数
% %            v(i,1)=0.5*v_limit(1);
%         end
%         if v(i,2)>v_limit(2)
%             v(i,2)=v_limit(2);
% %            v(i,2)=0.5*v_limit(2);
%         end
%         if v(i,3)>v_limit(3)
%             v(i,3)=v_limit(3);
% %             v(i,3)=0.5*v_limit(3);
%         end
%         if v(i,4)>v_limit(4)
%             v(i,4)=v_limit(4);
% %             v(i,4)=0.5*v_limit(4);
%         end
%         if v(i,5)>v_limit(5)
%             v(i,5)=v_limit(5);
% %             v(i,5)=0.5*v_limit(5);
%         end
%         if v(i,6)>v_limit(6)
%             v(i,6)=v_limit(6);
% %              v(i,6)=0.5*v_limit(6);
%         end
%         if v(i,7)>v_limit(7)
%             v(i,7)=v_limit(7);
% %             v(i,7)=0.5*v_limit(7);
%         end
%         if v(i,8)>v_limit(8)
%             v(i,8)=v_limit(8);
% %             v(i,8)=0.5*v_limit(8);
%         end
%         if v(i,9)>v_limit(9)
%             v(i,9)=v_limit(9);
% %             v(i,9)=0.5*v_limit(9);
%         end
%         if v(i,10)>v_limit(10)
%             v(i,10)=v_limit(10);
% %             v(i,10)=0.5*v_limit(10);
%         end
%         
%         x(i,:)=x(i,:)+v(i,:);
%         
% %         if(x(i,1)>xMax(1)||x(i,2)>xMax(2)||x(i,3)>xMax(3)...
% %                 ||x(i,4)>xMax(4)||x(i,5)>xMax(5)||x(i,6)>xMax(6)...
% %                 ||x(i,7)>xMax(7)||x(i,8)>xMax(8)||x(i,9)>xMax(9)||x(i,10)>xMax(10)...
% %                 ||x(i,1)<xMin(1)||x(i,2)<xMin(2)||x(i,3)<xMin(3)...
% %                 ||x(i,4)<xMin(4)||x(i,5)<xMin(5)||x(i,6)<xMin(6)...
% %                 ||x(i,7)<xMin(7)||x(i,8)<xMin(8)||x(i,9)<xMin(9)||x(i,10)<xMin(10)...
% %                 ||x(i,1)<=x(i,2))
% %             continue;
% %         end
% 
%         for index=1:10
%             if(x(i,index)>xMax(index))
%                 x(i,index)=xMax(index);
%             end
%             if(x(i,index)<xMin(index))
%                 x(i,index)=xMin(index);
%             end
%         end
% 
%         value=fitness(x(i,:));
%         
%         if ~isnan(value)    %判断value是不是Nan，如果是Nan，直接跳过
%             if value<p(i)   
%                 p(i)=value;  
%                 y(i,:)=x(i,:);
%                 if p(i)<local_opt        %%%%一定是p(i)而不是value
%                     pg=y(i,:);
%                     local_opt=p(i);
%                 end
%             end
%         end
%     end
%     disp("current: "+t+" total: "+M +" time= "+toc+" opt= "+local_opt);
%     result(t)=local_opt;
%     xlabel=1:length(result);
%     plot(xlabel,result);
%     pause(0.5);
%     if local_opt<0.01
%         break;
%     end
% end
% %%
% %结束
% save fv local_opt -ascii;
% xm=pg';
% fv=local_opt;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%function[xm,fv]=PSO(fitness,40,2,2,0.5,5000,9)
% % function[xm,fv]=PSO(fitness,N,c1,c2,w,M,D)
% % %待优化得目标函数fitness；
% % %%粒子数目N
% % %%学习因子C1\C2
% % %%惯性权重w
% % %%最大迭代次数M
% % %%问题得维数D
% % %%目标函数取最小值时得自变量值xm
% % %%目标函数最小值fv
% % format long;
% %  tic;
% % format long;
% % % w_ini=0.9;
% % % w_end=0.4;
% % 
% % xMax1=10^-3;   %%%%%%%%%Kr范围，单位m/s%%%%%%%%%
% % xMin1=10^-9;
% % x1=rand(N,1)*(xMax1-xMin1)+xMin1;
% % 
% % xMax2=10^-3;    %%%%%%%%%Kz范围,单位m/s%%%%%%%%%
% % xMin2=10^-9;
% % x2=rand(N,1)*(xMax2-xMin2)+xMin2;
% % %Kz应该小于0.1Kr
% % 
% % xMax3=0.01;     %%%%%%%%%S范围%%%%%%%%%
% % xMin3=0.000000001;
% % x3=rand(N,1)*(xMax3-xMin3)+xMin3;
% % 
% % xMax4=3;        %%%%%%%%%P范围%%%%%%%%%
% % xMin4=1;
% % x4=rand(N,1)*(xMax4-xMin4)+xMin4;
% % 
% % 
% % xMax5=1;     %%%%%%%%%C范围%%%%%%%%%
% % xMin5=1e-12;
% % x5=rand(N,1)*(xMax5-xMin5)+xMin5;
% % 
% % 
% % xMax6=1;     %%%%%%%%%alpha1范围%%%%%%%%%
% % xMin6=-5;
% % x6=rand(N,1)*(xMax6-xMin6)+xMin6;
% % 
% % xMax7=1;     %%%%%%%%%alpha2范围%%%%%%%%%
% % xMin7=-5;
% % x7=rand(N,1)*(xMax7-xMin7)+xMin7;
% % 
% % xMax8=1;     %%%%%%%%%alpha3范围%%%%%%%%%
% % xMin8=-5;
% % x8=rand(N,1)*(xMax8-xMin8)+xMin8;
% % 
% % xMax9=1;     %%%%%%%%%alpha4范围%%%%%%%%%
% % xMin9=-5;
% % x9=rand(N,1)*(xMax9-xMin9)+xMin9;
% % 
% % xMax10=1;     %%%%%%%%%alpha5范围%%%%%%%%%
% % xMin10=-5;
% % x10=rand(N,1)*(xMax10-xMin10)+xMin10;
% % 
% % xMax=[xMax1 xMax2 xMax3 xMax4 xMax5 xMax6 xMax7 xMax8 xMax9 xMax10];
% % xMin=[xMin1 xMin2 xMin3 xMin4 xMin5 xMin6 xMin7 xMin8 xMin9 xMin10];
% % x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10];
% % % xMax1=2;
% % % xMin1=0;
% % % x1=rand(N,1)*(xMax1-xMin1)+xMin1;
% % % xMax2=20;
% % % xMin2=0;
% % % x2=rand(N,1)*(xMax2-xMin2)+xMin2;
% % % xMax3=10;
% % % xMin3=0;
% % % x3=rand(N,1)*(xMax3-xMin3)+xMin3;
% % % xMax4=0.1;
% % % xMin4=0;
% % % x4=rand(N,1)*(xMax4-xMin4)+xMin4;
% % % xMax5=0.1;
% % % xMin5=0;
% % % x5=rand(N,1)*(xMax5-xMin5)+xMin5;
% % % xMax6=0.01;
% % % xMin6=0.000001;
% % % x6=rand(N,1)*(xMax6-xMin6)+xMin6;
% % % xMax7=0.01;
% % % xMin7=0.000001;
% % % x7=rand(N,1)*(xMax7-xMin7)+xMin7;
% % % xMax8=0.01;
% % % xMin8=0.000001;
% % % x8=rand(N,1)*(xMax8-xMin8)+xMin8;
% % % xMax9=0.01;
% % % xMin9=0.000001;
% % % x9=rand(N,1)*(xMax9-xMin9)+xMin9;
% % % x=[x1,x2,x3,x4,x5,x6,x7,x8,x9];
% % vMax1=5;
% % vMin1=-5;
% % v=rand(N,D)*(vMax1-vMin1)+xMin1;
% % d=1:N;
% % d=d';
% % dd=d;
% % % xMax=[xMax1 xMax2 xMax3 xMax4 xMax5 xMax6 xMax7 xMax8 xMax9];
% % [d,xMax]=ndgrid(d,xMax);
% % % xMin=[xMin1 xMin2 xMin3 xMin4 xMin5 xMin6 xMin7 xMin8 xMin9];
% % [dd,xMin]=ndgrid(dd,xMin);
% % 
% % p=zeros(N,1);
% % y=x;
% % pg=x(1,:);
% % local_opt=inf;
% % for i=1:N
% %     value=fitness(x(i,:));
% %     if isnan(value)
% %         p(i)=Inf;%p是处理后的结果矩阵
% %     else
% %     p(i)=value;
% %     y(i,:)=x(i,:);
% %         if value<local_opt    %local_opt 局部最优
% %             pg=x(i,:);
% %             local_opt=value;
% %         end
% %     end
% % end
% % 
% % 
% % % for i=1:N
% % %     value=fitness(x(i,:));
% % %     if isnan(value)==1
% % %         p(i)=Inf;
% % %     else
% % %     p(i)=value;
% % %     end
% % %     y(i,:)=x(i,:);
% % % end
% % % pg=x(N,:);
% % % %pg=x(N,:)
% % % for i=1:(N-1)
% % %     %for i=1:(N-1)
% % %     value1=fitness(x(i,:));
% % %     if isnan(value1)==1
% % %         comp1=Inf;
% % %     else comp1=value1;
% % %     end
% % %         if comp1<fitness(pg)
% % %             pg=x(i,:);
% % %         end
% % % end
% % result=zeros(1,M);
% % for t=1:M
% %     %for t=1:M
% %     for i=1:N
% %         %for i=1:N
% %         v(i,:)=w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
% %         x(i,:)=x(i,:)+v(i,:);
% %         
% %         bri1=x(i,:)-xMax(i,:);
% %         loc1=find(bri1>0);
% %         x(i,loc1)=0.5*xMax(i,loc1);
% %         
% %         bri2=x(i,:)-xMin(i,:);
% %         loc2=find(bri2<0);
% %         x(i,loc2)=xMin(i,loc2);
% %         
% %         value2=fitness(x(i,:));
% %     if isnan(value2)==1
% %         comp2=Inf;
% %     else comp2=value2;
% %     end
% %         if comp2<p(i)
% %             p(i)=comp2;
% %             y(i,:)=x(i,:);
% %         end
% %         if p(i)<local_opt
% %             pg=y(i,:);
% %             local_opt=p(i);
% %         end
% %     end
% %     disp("current: "+t+" total: "+M +" time= "+toc+" opt= "+local_opt);
% %     result(t)=local_opt;
% %     xlabel=1:length(result);
% %     plot(xlabel,result);
% %     pause(0.5);
% % %     pbest(t)=fitness(pg);
% %     if local_opt<0.01
% %         break;
% %     end
% % end
% % % save fv pbest -ascii;
% % xm=pg';
% % fv=local_opt;
%           
%         
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
