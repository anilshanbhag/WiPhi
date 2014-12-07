function [x,val,z]=optimize_iterative_2(peaks,X,lambda,z,exhaustive)
n_peaks=length(peaks);
L=repmat(3*pi./lambda',1,n_peaks);
D=repmat(peaks,length(lambda),1);
A=exp(-1j.*L.*D);
pA=pinv(A);
converged=false;
rng('shuffle');    
max_iter=4000;
iter=1;
if(exhaustive)
    v=zeros(2048,1);
    for i=0:2047
        p=dec2bin(i,11);
        z=(int32(p)*2-97)';
        z=(z>0)*2-1;
        x1=pA*(sqrt(X).*z);
        v(i+1)=norm((A*x1)-sqrt(X).*z);
    end
    [val,idx]=min(v);
    p=dec2bin(idx-1,11);
    z=(int32(p)*2-97)';
    z=(z>0)*2-1;
    x=pA*(sqrt(X).*z);
    fprintf(1,'\n Done with exhaustive');

else
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
            %fprintf(2,num2str(z'));
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
end

