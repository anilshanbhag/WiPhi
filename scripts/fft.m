d_vals=0:0.01:7.5;
for i=1:length(d_vals)
    DP(i)=sum(h.*exp(-1j*4*pi./lambda*d_vals(i)));   %lambda is a vector of frequencies and h is a vector of corresponding channels.
end
figure; plot(abs(DP));
