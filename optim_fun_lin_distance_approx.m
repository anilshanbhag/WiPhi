function u=optim_fun_lin_distance_approx(x,X,lambda,options)  


%X has to be a column vectors (square root is done inside the function)
%lambda has to a row vector

x=x.';
n_peaks=length(x);
p_factor=options.p_factor/2;

l1=ceil(length(lambda)/2);
i_s=0:1:2^(l1-1)-1;
p_s=dec2bin(i_s,l1-1);
z_s=(int32(p_s)*2-97)';
z_s=(z_s>0)*2-1;
z_s=[ones(1,length(i_s));z_s];
X_s=repmat(sqrt(X(1:l1)),1,length(i_s)).*z_s;
L=repmat(-1j*p_factor*pi./lambda',1,n_peaks);
[~,~,idx]=optimize_iterative_3(x,L(1:l1,:),X_s);    
z1=z_s(:,idx);

l2=length(lambda)-l1+1;
i_s=0:1:2^(l2-1)-1;
p_s=dec2bin(i_s,l2-1);
z_s=(int32(p_s)*2-97)';
z_s=(z_s>0)*2-1;
z_s=[ones(1,length(i_s));z_s];
X_s=repmat(sqrt(X([1,l1+1:length(lambda)])),1,length(i_s)).*z_s;
L=repmat(-1j*p_factor*pi./lambda',1,n_peaks);
[~,~,idx1]=optimize_iterative_3(x,L([1,l1+1:length(lambda)],:),X_s);    
z2=z_s(:,idx1);

z=[z1;z2(2:end)];

X_f=sqrt(X).*z;
[~,u]=optimize_iterative_3(x,L,X_f);


% alpha=options.alpha;
% delta=options.delta;
% l=length(lambda);
% x=x.';
% X=sqrt(X);
% L=repmat(2*pi*1./lambda',1,length(x));
% D=repmat(x,length(lambda),1);
% A1=[diag(X),exp(-1j.*L.*D)];
% h1= @(x1) optim_fun_lin(x1,A1,length(lambda),alpha);
% opt=optimoptions('fmincon','TolFun',1e-15,'TolCon',1e-15,'TolX',1e-15,'Display','off');
% x_opt=[randi(2,11,1)*2-3;rand(length(x),1)*5];        
% v=fmincon(h1,x_opt,[],[],[],[],[-1*ones(11,1);-10*ones(length(x),1)],[ones(11,1);10*ones(length(x),1)],[],opt); 
% z=sign(v(1:l));
% a=-pinv(A1(:,l+1:end))*(X.*z);
% u=h1([z;a]);
% 
% idx=randi(length(x));
% x_new=x;
% x_new(idx)=x(idx)+(randi(2)*2-3)*delta;
% L=repmat(2*pi*1./lambda',1,length(x_new));
% D=repmat(x_new,length(lambda),1);
% A1=[diag(X),exp(-1j.*L.*D)];
% h1= @(x1) optim_fun_lin(x1,A1,length(lambda),alpha);
% opt=optimoptions('fmincon','TolFun',1e-15,'TolCon',1e-15,'TolX',1e-15,'Display','off');
% x_opt=[randi(2,11,1)*2-3;rand(length(x),1)*5];        
% v=fmincon(h1,x_opt,[],[],[],[],[-1*ones(11,1);-10*ones(length(x_new),1)],[ones(11,1);10*ones(length(x_new),1)],[],opt); 
% z=sign(v(1:l));
% a=-pinv(A1(:,l+1:end))*(X.*z);
% u1=h1([z;a]);
% 
% 
% L=repmat(2*pi*1./lambda',1,length(x));
% D=repmat(x,length(lambda),1);
% A1=[diag(X),exp(-1j.*L.*D)];
% h1= @(x1) optim_fun_lin(x1,A1,length(lambda),alpha);
% a=-pinv(A1(:,l+1:end))*(X.*z);
% u2=h1([z;a]);
% 
% if(u2<u)
%     u=u2;
% end


end