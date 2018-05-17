function sValues=steps2goal(S,T)
sValues=[S,zeros(size(S,1),1)];
step=0;
sGoalIndex=112;
Tmod=T-eye(225);
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