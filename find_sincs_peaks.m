function [A,I]=find_sincs_peaks(X,lambda,options)

%X has to be a column vectors (square root is done inside the function)
%lambda has to a row vector

n_peaks=options.n_peaks;
res=options.res;
min_val=options.min_val;
max_val=options.max_val;
p_factor=options.p_factor/2;
n_vals=round((max_val-min_val)/res)+1;
A=zeros(n_vals^n_peaks,1);
I=zeros(n_vals^n_peaks,1);
i_s=0:1:2^(length(lambda)-1)-1;
p_s=dec2bin(i_s,length(lambda)-1);
z_s=(int32(p_s)*2-97)';
z_s=(z_s>0)*2-1;
z_s=[ones(1,length(i_s));z_s];
X_s=repmat(sqrt(X),1,length(i_s)).*z_s;
%X_s=repmat((X),1,length(i_s)).*z_s;
L=repmat(-1j*p_factor*pi./lambda',1,n_peaks);
x=zeros(1,n_peaks);
for i=0:1:n_vals^n_peaks-1
    x=min_val+floor(mod(i,n_vals.^(1:n_peaks))./n_vals.^(0:n_peaks-1))*res;
    [~,A(i+1),I(i+1)]=optimize_iterative_3(x,L,X_s);
    if(mod(i+1,100)==0)
        disp([i+1, n_vals^n_peaks]);
    end
end
if(n_peaks==2)
    A=reshape(A,n_vals,n_vals);
    I=reshape(I,n_vals,n_vals);
    x_vals=min_val:res:max_val;
    h=figure; hold on; surf(x_vals,x_vals,1./A,'EdgeColor','none');  
    colormap jet
elseif(n_peaks==1)
    x_vals=min_val:res:max_val;
    figure; plot(x_vals,1./A);
elseif(n_peaks==3)
    figure;
    A=reshape(A,[n_vals,n_vals,n_vals]);
    x_vals=min_val:res:max_val;
    for i=1:1:n_vals
        clf;
        hold on; surf(x_vals,x_vals,1./squeeze(A(:,:,i)),'EdgeColor','none');    
        pause(0.1);
    end
end



end

