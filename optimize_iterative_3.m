function [x,v,idx,res]=optimize_iterative_3(peaks,L,X_s)

D=repmat(peaks,size(L,1),1);
A=exp(L.*D);
pA=pinv(A);
x1=pA*(X_s);
v=sum(abs(A*x1-X_s).^2,1);
[v,idx]=min(v);
x=x1(:,idx);
res=X_s(:,idx)-A*x;    
end

