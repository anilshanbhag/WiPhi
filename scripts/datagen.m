freq=[5.5:0.02:5.7]*1e9; % Frequencies where you want channels
lambda=3e8./freq;
d=[2,4];                % Distances where the signals are coming from
hr=zeros(1,length(lambda));
h_f=zeros(length(d),length(lambda));   
 
for i=1:length(d)
    for j=1:length(lambda)       
        h_f(i,j)=exp(-1j*d(i)*2*pi/lambda(j));       
    end   
end
h_r=sum(h_f,1);
h=h_r.^2;         % We are measuring h^2
h=awgn(h,30);     % Add some Gaussian noise
