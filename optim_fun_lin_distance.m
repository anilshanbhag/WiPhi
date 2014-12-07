function u=optim_fun_lin_distance(x,X,lambda,options)       
%X has to be a column vectors (square root is done inside the function)
%lambda has to a row vector
x=x.';
n_peaks=length(x);
p_factor=options.p_factor/2;
i_s=0:1:2^(length(lambda)-1)-1;
p_s=dec2bin(i_s,length(lambda)-1);
z_s=(int32(p_s)*2-97)';
z_s=(z_s>0)*2-1;
z_s=[ones(1,length(i_s));z_s];
X_s=repmat(sqrt(X),1,length(i_s)).*z_s;
L=repmat(-1j*p_factor*pi./lambda',1,n_peaks);
[~,u]=optimize_iterative_3(x,L,X_s);    



end