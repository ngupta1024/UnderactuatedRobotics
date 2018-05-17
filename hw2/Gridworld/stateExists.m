function [sIndex, sBool]=stateExists(stateCheck, stateSpace)
Stemp=stateSpace==stateCheck;
sIndex=find(Stemp(:,1).*Stemp(:,2));
if sIndex
    sBool=true;
else
    sBool=false;
end
end