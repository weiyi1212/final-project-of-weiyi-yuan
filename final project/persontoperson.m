% % % simple person to person model
% % % Q=fg(r,r0,r1)t;A dose of Q=1 corresponds to a 63% probability of infection
% % % f is a constant[1,400](for hour), g(r) is function of distance r, t is time(hour)
% % % g(r,r0,r1)=1/r when r0<=r<=r1 , g(r)=1 when r<r0, g(r)=0 when r>r1
% % % 
% % % 
% % % a square of L*W point with each point ony can stand 1 person
% % %                   *
% % %                   *
% % %                   L,i 
% % %                   *
% % %                   *  *  W,j *  *
% % % pl length (width) of one point , to calculate distance between 2 people (r)
% % % calculate at time t= dt, 2*dt, 3*dt....
% % % people could move to a point around them at time t, with probability of pm
% % % Do not consider infected people release the virus
% % % The virus inhaled during the simulation time is not destroyed, has been accumulated
clc
clear
dt=5.0/60;
f=20.0;r0=1.0;r1=5.0;
L=20;W=20;
Ninfe=1;Nunin=179;pl=0.5;pm=0.5;
qinfet=9999;
nums=floor(r1/pl);%The maximum distance between the infected area and the infected person
Ttotal=45.0/60;Tnum=floor(Ttotal/dt);
infect=zeros(Ninfe,2);%Record information about infected people. Position information: Column 1 is the row number and column 2 is the column number;ilocation,jlocation
uninfe=zeros(Nunin,4);%Recording information on infected persons; Position information: Column 1 is the row number and column 2 is the column numberilocation,jlocation?3列是迭代时段吸入的病毒量，4列是累积吸入量
dQi=3;aQi=4;%
people=cell(L,W,Tnum+1);%Record the serial number of people at each grid point at each time. Column 1 is infected and column 2 is audience
people(:,:,:)={(zeros(1,2))};
Qviral=zeros(L,W,Tnum+1);%The cumulative amount of virus at each grid point at each time
total_inf=zeros(Tnum+1,1);%
% initial
index=randperm(L*W);
index_infe=index(1:Ninfe);
index_unin=index(Ninfe+1:Ninfe+Nunin);
[i,j]=ind2sub([L,W],index_infe);%Infected number
infect(:,1)=i;
infect(:,2)=j;

[i,j]=ind2sub([L,W],index_unin);%Audience number
uninfe(:,1)=i;
uninfe(:,2)=j;

for i=1:Ninfe
    people{infect(i,1),infect(i,2),1}(1,1)=i;
end
for i=1:Nunin
    people{uninfe(i,1),uninfe(i,2),1}(1,2)=i;
end
% % % Virus quantity at point (iL,iW) in DT in 1 iteration period
% % % Because the number of infected people is low, so the circulation of infected people
uninfe(:,dQi)=0.0;
for iN=1:Ninfe
    ilocat=infect(iN,1);
    jlocat=infect(iN,2);
    iLs=max(ilocat-nums,1);iLe=min(ilocat+nums,L);%The maximum and minimum of the lattice sequence number
    iWs=max(jlocat-nums,1);iWe=min(jlocat+nums,W);
    for iL=iLs:iLe
    for iW=iWs:iWe%euler
        index_in=people{iL,iW,1}(1,2);%The audience the serial number       
        if index_in==0
            continue
        else
            r=pl*sqrt((iL-ilocat)^2+(iW-jlocat)^2);%distance between infected people and the uninfected
            dQ=f*dt*g(r,r0,r1);%time step dt Amount of virus inhaled           
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
% % %People move randomly
move=1:Ninfe+Nunin;%Serial numbers of all the people to be moved
for i=1:Ninfe+Nunin% 
    if ismember(i,move)&&randsrc(1,1,[0,1;1-pm,pm])%If the i th person wants to move
        if i<=Ninfe%if This man is infected
            iindex=i;
            ilocat=infect(iindex,1);
            jlocat=infect(iindex,2);
        else%audience?           
            iindex=i-Ninfe;
            ilocat=uninfe(iindex,1);
            jlocat=uninfe(iindex,2);
        end
        ilocat0=randsrc(1,1,[1,2,3,4;0.25,0.25,0.25,0.25]);%movement direction.up, down, left, right
        if ilocat0==1%up?           
            ilocat0=min(ilocat+1,L);jlocat0=jlocat;
        elseif ilocat0==2%right?           
            ilocat0=ilocat;jlocat0=min(jlocat+1,W);
        elseif ilocat0==3%down?           
            ilocat0=max(ilocat-1,1);jlocat0=jlocat;  
        else%left?           
            ilocat0=ilocat;jlocat0=max(jlocat-1,1); 
        end
        a=people{ilocat0,jlocat0,it-1};%  Update the sequence number of the grid point?
        people{ilocat0,jlocat0,it}=people{ilocat,jlocat,it-1};
        people{ilocat,jlocat,it}=a;
% % %   Update people's location information
        if i<=Ninfe%renew ith person
            infect(iindex,1)=ilocat0;
            infect(iindex,2)=jlocat0;
        else
            uninfe(iindex,1)=ilocat0;
            uninfe(iindex,2)=jlocat0;
        end
% % %   delete i in move
        move(move==i)=[];        
% % % If there's someone at the target, update that person's location?        
        if a(1,2)>0 %The target is the audience
            j=a(1,2);
            uninfe(j,1)=ilocat;
            uninfe(j,2)=jlocat;
            j=j+Ninfe;
            move(move==j)=[];% %  delete j in move
        elseif a(1,1)>0%The target is the infected        
            j=a(1,1);
            infect(j,1)=ilocat;
            infect(j,2)=jlocat;
            move(move==j)=[];%  delete j in move
        end
    else%if do not move        
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




