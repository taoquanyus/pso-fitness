function[xm,fv]=new_PSO(fitness,N,c1,c2,w_ini,w_end,M,D)
%% %λ������
tic;
xMax=[1e-2,1e-2,1e-2,3,1,3,3,3,3,3];
% xMax=[1e-2,1e-2,1e-2,3,1,3];
xMin=[1e-9,1e-9,1e-9,1,1e-14,-3,-3,-3,-3,-3];
% xMin=[1e-9,1e-9,1e-9,1,1e-14,-3];

%% %�ٶ�����
v_index=0.15;
vMax=v_index*xMax;
vMin=-vMax;

%% %��ʼ��λ��,�ٶ�

x_cur=rand(N,D).*(xMax-xMin)+xMin;
% x_cur(1,:)=[0.000585764066514 0.000123132275533 0.005415238111098 1.408037862239064 0.414684292768394 -2.43 -1.58 -1.58 -1.65 -1.71];
v=rand(N,D).*(vMax-vMin)+vMin;

%% %�״ε���
results=zeros(N,1);%�������
pBest=x_cur;
gBest=x_cur(1,:);
pBest_result=inf;%�ֲ����Ž��ʼ��
gBest_result=inf;%ȫ�����Ž��ʼ��

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

%% %M������
for times=1:M
    w=(w_ini-w_end)*(M-times)/M+w_end;
    vLarge_exceed=0;
    vSmall_exceed=0;
    xLarge_exceed=0;
    for row=1:N
        v(row,:)=w*v(row,:)+c1*rand*(pBest(row,:)-x_cur(row,:))+c2*rand*(gBest-x_cur(row,:));
        
        %�ٶ�����
        for temp=1:D
            if v(row,temp)>vMax(temp)
                v(row,temp)=vMax(temp);
                vLarge_exceed=vLarge_exceed+1;
            end
            if v(row,temp)<vMin(temp)
                v(row,temp)=vMin(temp);
                vSmall_exceed=vSmall_exceed+1;
            end
        end
        
        %λ������
        x_cur(row,:)=x_cur(row,:)+v(row,:);
        
        if x_cur(row,2)>x_cur(row,1)
            x_cur(row,2)=x_cur(row,1);
        end
        
        for temp=1:D
            if x_cur(row,temp)>xMax(temp)
                x_cur(row,temp)=xMax(temp);
                xLarge_exceed=xLarge_exceed+1;
            end
            if x_cur(row,temp)<xMax(temp)
                x_cur(row,temp)=xMin(temp);
                xSmall_exceed=xSmall_exceed+1;
            end
        end
        
        %����
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
    disp("current: "+times+" total: "+M +" time= "+toc+" opt= "+gBest_result+" V_exceed: "+vLarge_exceed+" "+vSmall_exceed+" x_exceed "+xLarge_exceed+" "+xSmall_exceed);
%     disp(x_cur');
    result(times)=gBest_result;
    plot(result,'o');
    pause(0.5);
    
    if gBest_result<0.1
        break;
    end
end

%%
%����
save fv gBest_result -ascii;
xm=gBest';
fv=gBest_result;