% % % simple person to person model
% % % Q=fg(r,r0,r1)t;Q=1�ļ�����Ӧ63%�ĸ�Ⱦ����
% % % f is a constant[1,400](for hour), g(r) is function of distance r, t is time(hour)
% % % g(r,r0,r1)=1/r2 when r0<=r<=r1 , g(r)=1 when r<r0, g(r)=0 when r>r1
% % % 
% % % 
% % % a square of L*W point with each point ֻ��վ?1 person
% % %                   *
% % %                   *
% % %                   L,i 
% % %                   *
% % %                   *  *  W,j *  *
% % % pl length (width) of one point , to calculate distance between 2 people (r)
% % % calculate at time t= dt, 2*dt, 3*dt....
% % % people could move to a point around them at time t, with probability of pm
% % % �����Ǹ�Ⱦ������ͷŲ���
% % % ģ��ʱ��������Ĳ�����������һֱ�ۻ�
clc
clear
dt=5.0/60;
f=20.0;r0=1.0;r1=5.0;
L=20;W=20;
Ninfe=1;Nunin=179;pl=0.5;pm=0.5;
qinfet=9999;
nums=floor(r1/pl);%��Ⱦ��Χ���Ⱦ�ߵ�����������?
Ttotal=45.0/60;Tnum=floor(Ttotal/dt);
infect=zeros(Ninfe,2);%��¼��Ⱦ�ߵ���Ϣ?λ����Ϣ��1��������ţ�2���������ilocation,jlocation
uninfe=zeros(Nunin,4);%��¼��Ⱦ�ߵ���Ϣ?λ����Ϣ��1��������ţ�2���������ilocation,jlocation?3���ǵ���ʱ������Ĳ�������4�����ۻ�������
dQi=3;aQi=4;%
people=cell(L,W,Tnum+1);%��¼ÿ��ʱ��ÿ����㴦���˵���ţ�1���Ǹ�Ⱦ�ߣ�2���ǹ���
people(:,:,:)={(zeros(1,2))};
Qviral=zeros(L,W,Tnum+1);%ÿ��ʱ��ÿ����㴦���ۻ�������
total_inf=zeros(Tnum+1,1);%
% initial
index=randperm(L*W);
index_infe=index(1:Ninfe);
index_unin=index(Ninfe+1:Ninfe+Nunin);
[i,j]=ind2sub([L,W],index_infe);%��Ⱦ�����
infect(:,1)=i;
infect(:,2)=j;

[i,j]=ind2sub([L,W],index_unin);%�������?
uninfe(:,1)=i;
uninfe(:,2)=j;

for i=1:Ninfe
    people{infect(i,1),infect(i,2),1}(1,1)=i;
end
for i=1:Nunin
    people{uninfe(i,1),uninfe(i,2),1}(1,2)=i;
end
% % % 1������ʱ��dt�� point (iL ,iW)���Ĳ�����?
% % % ��Ϊ��Ⱦ���٣�����ѭ����Ⱦ��
uninfe(:,dQi)=0.0;
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    iLs=max(ilocat-nums,1);iLe=min(ilocat+nums,L);%�����ŵļ���ֵ��Сֵ
    iWs=max(jlocat-nums,1);iWe=min(jlocat+nums,W);
    for iL=iLs:iLe
    for iW=iWs:iWe%ŷ��
        index_in=people{iL,iW,1}(1,2);%�������?        
        if index_in==0
            continue
        else
            r=pl*sqrt((iL-ilocat)^2+(iW-jlocat)^2);%distance between infected people and the uninfected
            dQ=f*dt*g(r,r0,r1);%time step dt �ڲ���������?            
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
% % %������ƶ�λ��?
move=1:Ninfe+Nunin;%���д��ƶ��˵����
for i=1:Ninfe+Nunin% 
    if ismember(i,move)&&randsrc(1,1,[0,1;1-pm,pm])%�����i����Ҫ�ƶ�
        if i<=Ninfe%������Ǹ�Ⱦ��
            iindex=i;
            ilocat=infect(iindex,1);
            jlocat=infect(iindex,2);
        else%����?            
            iindex=i-Ninfe;
            ilocat=uninfe(iindex,1);
            jlocat=uninfe(iindex,2);
        end
        ilocat0=randsrc(1,1,[1,2,3,4;0.25,0.25,0.25,0.25]);%�ƶ�������������
        if ilocat0==1%��?           
            ilocat0=min(ilocat+1,L);jlocat0=jlocat;
        elseif ilocat0==2%��?           
            ilocat0=ilocat;jlocat0=min(jlocat+1,W);
        elseif ilocat0==3%��?           
            ilocat0=max(ilocat-1,1);jlocat0=jlocat;  
        else%��?           
            ilocat0=ilocat;jlocat0=max(jlocat-1,1); 
        end
        a=people{ilocat0,jlocat0,it-1};%  ���¸�㴦�˵����?
        people{ilocat0,jlocat0,it}=people{ilocat,jlocat,it-1};
        people{ilocat,jlocat,it}=a;
% % %   �����˵�λ����Ϣ 
        if i<=Ninfe%renew ith person
            infect(iindex,1)=ilocat0;
            infect(iindex,2)=jlocat0;
        else
            uninfe(iindex,1)=ilocat0;
            uninfe(iindex,2)=jlocat0;
        end
% % %   delete i in move
        move(move==i)=[];        
% % % ���Ŀ�괦���ˣ���������˵�λ����Ϣ?        
        if a(1,2)>0%Ŀ�괦�ǹ���
            j=a(1,2);
            uninfe(j,1)=ilocat;
            uninfe(j,2)=jlocat;
            j=j+Ninfe;
            move(move==j)=[];% %  delete j in move
        elseif a(1,1)>0%Ŀ�괦�Ǹ�Ⱦ��?            
            j=a(1,1);
            infect(j,1)=ilocat;
            infect(j,2)=jlocat;
            move(move==j)=[];%  delete j in move
        end
    else%������ƶ�?        
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




