function [x,val]=optimize_iterative_4(peaks,X,lambda)
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
%z=[-1,1,-1,1,-1,-1,-1,-1,-1,1,-1]';
z=z.*r;    
max_iter=5000;
iter=1;
while ~converged
    X_new=sqrt(X).*z;
    x=pA*X_new;
    r=randi(length(z));
    found_lower=false;
    goal=norm(A*x-X_new);
    for i=1:1:length(z)
        idx=mod(i+r,length(z));
        if(idx==0)
            idx=length(z);
        end
        z_new=z;
        z_new(idx)=z_new(idx)*-1;
        v=norm(A*(pA*(sqrt(X).*z_new))-sqrt(X).*z_new);
        if(v<goal)
            found_lower=true;
            z=z_new;
            break
        end
        
    end
    
    if(~found_lower)
        val=norm(A*x-X_new);
        converged=true;
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

