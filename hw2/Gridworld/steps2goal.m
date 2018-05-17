function sValues=steps2goal(S,T)
sGoal=[4.5,4.5];
sValues=[S,zeros(size(S,1),1)];
step=0;
Stemp=S==sGoal;
sGoalIndex=find(Stemp(:,1).*Stemp(:,2));
Tmod=T-eye(100);
for i=1:size(S,1)
    sNewGoalIndex=[];
    sValues(sGoalIndex,3)=step;
%     text(sValues(sGoalIndex,1),sValues(sGoalIndex,2),num2str(step));
    for j=1:size(sGoalIndex,2)
        sNewGoalIndex=[sNewGoalIndex,find(Tmod(sGoalIndex(j),:))];
        Tmod(:,sGoalIndex(j))=zeros(size(S,1),1);
    end  
    sGoalIndex=unique(sNewGoalIndex);
    step=step+1;
    if Tmod==zeros(size(T))
        break;
    end
end
end