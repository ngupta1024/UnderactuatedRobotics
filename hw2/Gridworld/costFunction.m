function costFunction(N, stateSpace, actionSpace, f)

% setting up the environment
sGoalIndex=100;
obstacles=false;

%costs
phi=ones(size(stateSpace,1),1);
g=ones(size(stateSpace,1),1); %independent of actions
phi(sGoalIndex)=0;
g(sGoalIndex)=0;

if(obstacles)
    sObstacleIndex=[3 4 13 14 9 10 19 20 29 30 23 24 93 94 8 18 28 38 39 40 83 84 45 46 47 55 56 57 65 66 67];
    phi(sObstacleIndex)=14;
    g(sObstacleIndex)=14;
    scatter(stateSpace(sObstacleIndex,1),stateSpace(sObstacleIndex,2),500,'y','filled');
end

J=phi;
J_new=J;
policy=ones(size(stateSpace,1),1)*5;
actions={'D','U','L','R','S'};
for i=1:N
    for currStateIndex=1:size(stateSpace,1)
        labels(currStateIndex)=text(stateSpace(currStateIndex,1),stateSpace(currStateIndex,2),num2str(J_new(currStateIndex)));
        dirs(currStateIndex)=text(stateSpace(currStateIndex,1)-0.4,stateSpace(currStateIndex,2)+0,actions{policy(currStateIndex)});
        for currAction=1:size(actionSpace,2)
            nextState=f(currStateIndex,currAction);
            if nextState~=0
                J_update(currAction)=g(currStateIndex)+ J(nextState);
            else 
                J_update(currAction)=Inf;
            end
        end
        [J_new(currStateIndex),policy(currStateIndex)]=min(J_update);
    end
    if J==J_new
        break;
    else
        J=J_new;
        pause
        delete(labels(:));
        delete(dirs(:));
end
end

