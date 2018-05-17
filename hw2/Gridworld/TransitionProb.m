function T=TransitionProb(x,y)
%we have 5 actions for every state
%length(x)=number of states

%move up-y+1
%move down=y-1
%move left=x-1
%move right=x+1
%stay put
points=[x',y'];
%if an action leads to a state that is in the state space
%100x100-very sparse
T=eye(length(x),length(x));

for i=1:length(x)
    currState=[x(i),y(i)];
    newStates=[x(i)-1,y(i);x(i)+1,y(i);x(i),y(i)-1;x(i),y(i)+1];
    for k=1:size(newStates,1)
        bool=points==newStates(k,:);
        T(i,:)=T(i,:)+(bool(:,1).*bool(:,2))';      
    end
end
end