function [x,val]=optimize_iterative_1(peaks,X,lambda)
z=ones(size(X));
n_peaks=length(peaks);
L=repmat(3*pi./lambda',1,n_peaks);
D=repmat(peaks,length(lambda),1);
A=exp(-1j.*L.*D);
pA=pinv(A);
converged=false;
rng('shuffle');
r=randi(2,size(X))*2-3;
%01010000010
z=[-1,1,-1,1,-1,-1,-1,-1,-1,1,-1]';
%z=z.*r;    
max_iter=4000;
iter=1;
while ~converged
    X_new=sqrt(X).*z;
    x=pA*X_new;
    v=zeros(length(z),1);
    for i=1:1:length(z)
        z_new=z;
        z_new(i)=z_new(i)*-1;
        v(i)=norm(A*(pA*(sqrt(X).*z_new))-sqrt(X).*z_new);
        
    end
    [l,idx]=min(v);
    if(l<norm(A*x-X_new))
        z(idx)=-z(idx);
    else
        val=norm(A*x-X_new);
        converged=true;
        fprintf(2,num2str(z'));
    end
   % disp([iter norm((A*x).^2-X)]);
    if(iter>max_iter)
        converged=true;
        fprintf(2,'\n Stopping due to number of iterations');
    else
        iter=iter+1;
    end
end
    
end

