function [policy,V,S,A]=policyIteration

close all

% states are mxn dimensional vector
m=linspace(0,6.2832,101);
n=linspace(-2*pi,2*pi,51);
[x1,x2]=ndgrid(m,n);%mxn
x1=reshape(x1,[1,size(x1,1)*size(x1,2)]);
x2=reshape(x2,[1,size(x2,1)*size(x2,2)]);
S=[x1; x2]';

%actions
A=linspace(-5,5,101);

%initialize policy and value function
policy=zeros(size(S,1),1);
V=zeros(size(S,1),1);
theta=0.1;
gamma=0.9; %only for value to converge; not used in policy iteration

for IterPI=1:1000
    %policy evaluation
    for iterPE=1:100
        delta=0;
        for s=1:size(S,1)
            temp=V(s);
            V(s)=costFunction(S(s,:),policy(s))+gamma*ValueFunction(V,transitionFunc(s, find(A==policy(s)),S,A,n),S);
            delta=max(delta,norm(temp-V(s)));
        end
        if delta<theta
            break;
        end
    end
    
    %policy improvement
    policyStable=true;
    for s=1:size(S,1)
        temp=policy(s);
        for a=1:length(A)
            valueList(a)= costFunction(s,A(a))+1*ValueFunction(V,transitionFunc(s,a,S,A,n),S);
        end
        [~, idx]=min(valueList); %choose an action that minimizes the cost
        policy(s)=A(idx);
        if temp~=policy(s)
            policyStable=false;
        end
    end
    if policyStable==true
        break;
    end
end

Plmap=reshape(policy,length(m),length(n));
surf(Plmap')
xlim([1 101])
ylim([1 51])
set(gca,'XTick',[1:10:101])
set(gca,'XTickLabel',string(m(1:10:101)))
set(gca,'YTick',[1:10:51])
set(gca,'YTickLabel',string(n(1:10:51)))
view(2);
title('Policy Map');
colorbar;

figure;
Vmap=reshape(V,length(m),length(n));
surf(Vmap');
xlim([1 101])
ylim([1 51])
set(gca,'XTick',[1:10:101])
set(gca,'XTickLabel',string(m(1:10:101)))
set(gca,'YTick',[1:10:51])
set(gca,'YTickLabel',string(n(1:10:51)))
view(2);
title('Value Function Map')
colorbar;

        
%% Functions
    function VInterp=ValueFunction(V,Snext,S)
        if any(S(:,1)==Snext(1)) && any(S(:,2)==Snext(2))
            VInterp=V(find(S(:,1)==Snext(1) & S(:,2)==Snext(2)));
        else
            %find Si1 and Si2 wrt Snext(1) and Snext(2)
            greatervals=find(m>Snext(1));
            Si1=greatervals(1)-1;
            greatervals=find(n>Snext(2));
            Si2=greatervals(1)-1;
            Si1Si2=find(S(:,1)==m(Si1) & S(:,2)==n(Si2));
            Si1plus1Si2=find(S(:,1)==m(Si1+1) & S(:,2)==n(Si2));
            Si1Si2plus1=find(S(:,1)==m(Si1) & S(:,2)==n(Si2+1));
            Si1plus1Si2plus1=find(S(:,1)==m(Si1+1) & S(:,2)==n(Si2+1));
            
            %find alpha1 and 2
            alpha1=1-(Snext(1)-S(Si1Si2,1))/(m(Si1+1)-m(Si1));
            alpha2=1-(Snext(2)-S(Si1Si2,2))/(n(Si2+1)-n(Si2));
            
            VInterp=alpha1*alpha2*V(Si1Si2)+(alpha1)*(1-alpha2)*V(Si1Si2plus1)+...
                (1-alpha1)*alpha2*V(Si1plus1Si2)+(1-alpha1)*(1-alpha2)*V(Si1plus1Si2plus1);
        end
    end

    function J=costFunction(x,u)
        xd=[pi;0];
        Q=eye(2);
        R=1;
        J=(x'-xd)'*Q*(x'-xd)+u'*R*u;
    end

    function sNext=transitionFunc(s,a,S,A,n)
        %use simple euler integration
        %y1=y0+h*f(y0)
        dx1dt=S(s,2);
        dx2dt=-S(s,2)-sin(S(s,1))+A(a);
        sNext=zeros(1,2);
        sNext(1)=S(s,1)+0.01*dx1dt;
        if abs(S(s,2))==max(n)
            sNext(2)=min(n);
        else
            sNext(2)=S(s,2)+0.01*dx2dt;
        end
        sNext(1)=wrapTo2Pi(sNext(1));  
    end
end