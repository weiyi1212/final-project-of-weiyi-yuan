% % % simple person to person model
% % % Q=fg(r,r0,r1)t;Q=1的剂量对应63%的感染概率
% % % f is a constant[1,400](for hour), g(r) is function of distance r, t is time(hour)
% % % g(r,r0,r1)=1/r2 when r0<=r<=r1 , g(r)=1 when r<r0, g(r)=0 when r>r1
% % % 
% % % 
% % % a square of L*W point with each point 只能站?1 person
% % %                   *
% % %                   *
% % %                   L,i 
% % %                   *
% % %                   *  *  W,j *  *
% % % pl length (width) of one point , to calculate distance between 2 people (r)
% % % calculate at time t= dt, 2*dt, 3*dt....
% % % people could move to a point around them at time t, with probability of pm
% % % 不考虑感染后的人释放病毒
% % % 模拟时间内吸入的病毒不被消灭，一直累积
clc
clear
dt=5.0/60;
f=20.0;r0=1.0;r1=5.0;
L=20;W=20;
Ninfe=1;Nunin=179;pl=0.5;pm=0.5;
qinfet=9999;
nums=floor(r1/pl);%感染范围距感染者的最大距离格点数?
Ttotal=45.0/60;Tnum=floor(Ttotal/dt);
infect=zeros(Ninfe,2);%记录感染者的信息?位置信息：1列是行序号，2列是列序号ilocation,jlocation
uninfe=zeros(Nunin,4);%记录感染者的信息?位置信息：1列是行序号，2列是列序号ilocation,jlocation?3列是迭代时段吸入的病毒量，4列是累积吸入量
dQi=3;aQi=4;%
people=cell(L,W,Tnum+1);%记录每个时刻每个格点处的人的序号，1列是感染者，2列是观众
people(:,:,:)={(zeros(1,2))};
Qviral=zeros(L,W,Tnum+1);%每个时刻每个格点处的累积病毒量
total_inf=zeros(Tnum+1,1);%
% initial
index=randperm(L*W);
index_infe=index(1:Ninfe);
index_unin=index(Ninfe+1:Ninfe+Nunin);
[i,j]=ind2sub([L,W],index_infe);%感染者序号
infect(:,1)=i;
infect(:,2)=j;

[i,j]=ind2sub([L,W],index_unin);%观众序号?
uninfe(:,1)=i;
uninfe(:,2)=j;

for i=1:Ninfe
    people{infect(i,1),infect(i,2),1}(1,1)=i;
end
for i=1:Nunin
    people{uninfe(i,1),uninfe(i,2),1}(1,2)=i;
end
% % % 1个迭代时段dt内 point (iL ,iW)处的病毒量?
% % % 因为感染者少，所以循环感染者
uninfe(:,dQi)=0.0;
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    iLs=max(ilocat-nums,1);iLe=min(ilocat+nums,L);%格点序号的极大值极小值
    iWs=max(jlocat-nums,1);iWe=min(jlocat+nums,W);
    for iL=iLs:iLe
    for iW=iWs:iWe%欧拉
        index_in=people{iL,iW,1}(1,2);%观众序号?        
        if index_in==0
            continue
        else
            r=pl*sqrt((iL-ilocat)^2+(iW-jlocat)^2);%distance between infected people and the uninfected
            dQ=f*dt*g(r,r0,r1);%time step dt 内病毒吸入量?            
            uninfe(index_in,dQi)=uninfe(index_in,dQi)+dQ;
            uninfe(index_in,aQi)=uninfe(index_in,aQi)+dQ;
        end
    end
    end
end
Qviral(:,:,1)=0.0;
for iN=1:Nunin
    ilocat=uninfe(iN,1);
    jlocat=uninfe(iN,2);
    Qviral(ilocat,jlocat,1)=uninfe(iN,aQi);
end
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    Qviral(ilocat,jlocat,1)=qinfet;
end
total_inf(1,1)=length(find(Qviral(:,:,1)*0.63>1.0));
for it=2:Tnum+1
people(:,:,it)=people(:,:,it-1);
% % %人随机移动位置?
move=1:Ninfe+Nunin;%所有待移动人的序号
for i=1:Ninfe+Nunin% 
    if ismember(i,move)&&randsrc(1,1,[0,1;1-pm,pm])%如果第i个人要移动
        if i<=Ninfe%这个人是感染者
            iindex=i;
            ilocat=infect(iindex,1);
            jlocat=infect(iindex,2);
        else%观众?            
            iindex=i-Ninfe;
            ilocat=uninfe(iindex,1);
            jlocat=uninfe(iindex,2);
        end
        ilocat0=randsrc(1,1,[1,2,3,4;0.25,0.25,0.25,0.25]);%移动方向，上下左右
        if ilocat0==1%上?           
            ilocat0=min(ilocat+1,L);jlocat0=jlocat;
        elseif ilocat0==2%右?           
            ilocat0=ilocat;jlocat0=min(jlocat+1,W);
        elseif ilocat0==3%下?           
            ilocat0=max(ilocat-1,1);jlocat0=jlocat;  
        else%左?           
            ilocat0=ilocat;jlocat0=max(jlocat-1,1); 
        end
        a=people{ilocat0,jlocat0,it-1};%  更新格点处人的序号?
        people{ilocat0,jlocat0,it}=people{ilocat,jlocat,it-1};
        people{ilocat,jlocat,it}=a;
% % %   更新人的位置信息 
        if i<=Ninfe%renew ith person
            infect(iindex,1)=ilocat0;
            infect(iindex,2)=jlocat0;
        else
            uninfe(iindex,1)=ilocat0;
            uninfe(iindex,2)=jlocat0;
        end
% % %   delete i in move
        move(move==i)=[];        
% % % 如果目标处有人，更新这个人的位置信息?        
        if a(1,2)>0%目标处是观众
            j=a(1,2);
            uninfe(j,1)=ilocat;
            uninfe(j,2)=jlocat;
            j=j+Ninfe;
            move(move==j)=[];% %  delete j in move
        elseif a(1,1)>0%目标处是感染者?            
            j=a(1,1);
            infect(j,1)=ilocat;
            infect(j,2)=jlocat;
            move(move==j)=[];%  delete j in move
        end
    else%如果不移动?        
        continue
    end
end

uninfe(:,dQi)=0.0;%
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    iLs=max(ilocat-nums,1);iLe=min(ilocat+nums,L);
    iWs=max(jlocat-nums,1);iWe=min(jlocat+nums,W);
    for iL=iLs:iLe
    for iW=iWs:iWe
        index_in=people{iL,iW,it}(1,2);%
        if index_in==0
            continue
        else
            r=pl*sqrt((iL-ilocat)^2+(iW-jlocat)^2);
            dQ=f*dt*g(r,r0,r1);%
            uninfe(index_in,dQi)=uninfe(index_in,dQi)+dQ;
            uninfe(index_in,aQi)=uninfe(index_in,aQi)+dQ;
        end
    end
    end
end
Qviral(:,:,it)=0.0;
for iN=1:Nunin
    ilocat=uninfe(iN,1);
    jlocat=uninfe(iN,2);
    Qviral(ilocat,jlocat,it)=uninfe(iN,aQi);
end
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    Qviral(ilocat,jlocat,it)=qinfet;
end
total_inf(it,1)=length(find(Qviral(:,:,it)>1.0));
end




