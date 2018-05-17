function learning
actionSpace=[1,2,3,4,5];%down=1;up=2;left=3;right=4;stayput=5;
[x,y]=CreateGrid;
T=TransitionProb(x,y);
S=[x',y'];
stateSpace=steps2goal(S,T); %with 3rd column with the value
f=stateTransition(stateSpace, actionSpace);
costFunction(100, stateSpace, actionSpace, f);
end  