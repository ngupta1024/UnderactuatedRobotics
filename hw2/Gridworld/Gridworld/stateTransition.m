function f=stateTransition(stateSpace, actionSpace)

%f is a matrix with states,s and actions,a-gives out s'
f=zeros(size(stateSpace,1),length(actionSpace));
for currStateIndex=1:size(stateSpace,1)
    f(currStateIndex,5)=currStateIndex;%stayput
    for currAction=1:length(actionSpace)-1
        if currAction==1
            sTransition=stateSpace(currStateIndex,1:2)+[0,-1];   
        elseif currAction==2
            sTransition=stateSpace(currStateIndex,1:2)+[0,1];
        elseif currAction==3
            sTransition=stateSpace(currStateIndex,1:2)+[-1,0];
        elseif currAction==4
            sTransition=stateSpace(currStateIndex,1:2)+[1,0];
        end
        [sTransIndex, existBool]=stateExists(sTransition, stateSpace(:,1:2));
        if existBool
                f(currStateIndex,currAction)=sTransIndex;
        end
    end
end
end
