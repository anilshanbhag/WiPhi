clear all
close all
d_vals=0:0.001:7.5; %Distance values to consider for distance profiles


% You can use either 24 channels or 11 channels.11 channels are regularly
% spaced and hence give a better distance profile. 
freq=(5.5:0.02:5.7)*1e9;
chan=[100:4:140]+1;

%freq=[5180,5200,5220,5240,5260,5280,5300,5320,5500,5520,5540,5560,5580,5600,5620,5640,5660,5680,5700,5745,5765,5785,5805,5825]*1e6;
%chan=[36:4:64,100:4:140,149:4:165]+1;

lambda=3e8./freq;
lambda_vals=1./lambda;
D=repmat(d_vals',1,length(lambda_vals));
L=repmat(lambda_vals,length(d_vals),1);
p_factor=4;
M=exp(1j.*p_factor*pi.*L.*D);
files_all=[1:15,18,19,16,17,20,21,24,25,22,23,26:37];
ground_truth=[18,34,51,72,80,87,92,97,104,112,112,119,119,126,126,134,134,139,139,142,142,147,147,150,150,155,155,158,158,166,166,167,167,176,176,185,185]*0.0254+1.83;

for file_idx=1:length(files_all)
    filename{1}=sprintf('Data\\3Dec_m1_%d.dat',files_all(file_idx));
    filename{2}=sprintf('Data\\3Dec_m2_%d.dat',files_all(file_idx));
    ant_num=1;
    i=1;[channels{i}, timestamp{i},packet_ids{i},MACs{i},sec{i},usec{i}]=process_trace_channel(filename{i},[]);
    i=2;[channels{i}, timestamp{i},packet_ids{i},MACs{i},sec{i},usec{i}]=process_trace_channel(filename{i},[]);               
    
    %% Channel Computation: Pre processing step. You can skip this.
    for ch=1:1:length(chan)
        if(chan(ch)<150)
            sub = [-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,28];
        else
            sub= [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58];
        end

        idx1=find(MACs{1}==chan(ch));
        idx2=find(MACs{2}==chan(ch));
        if(length(idx1)~=length(idx2))
            j=1;
            valid=false(size(idx1));
            p1=packet_ids{1}(idx1);
            p2=packet_ids{2}(idx2);
            for i=1:1:length(idx2)
                while(j<=length(p1) && p1(j)~=p2(i) && j<=length(idx1))
                    j=j+1;
                end
                if(j>length(p1))
                    idx2=idx2(1:end-1);
                    break;
                end
                valid(j)=true;
            end
            idx1=idx1(valid);
        end
        ant_second=false;
        for i=1:1:length(idx1)
            if(~ant_second)
                m=squeeze(angle(channels{1}(ant_num,idx1(i),:)));
                m1=squeeze(angle(channels{2}(1,idx2(i),:)));     
            else
                m=squeeze(angle(channels{1}(2,idx1(i),:)));
                m1=squeeze(angle(channels{2}(2,idx2(i),:)));            
            end
            slope=regress(unwrap(m),[sub', ones(length(sub),1)]);
            if(~ant_second)
                h1=squeeze(channels{1}(ant_num,idx1(i),:)).*exp(-1j*slope(1)*sub');
            else
                h1=squeeze(channels{1}(2,idx1(i),:)).*exp(-1j*slope(1)*sub');
            end
            slope=regress(unwrap(m1),[sub', ones(length(sub),1)]);
            if(~ant_second)
                h2=squeeze(channels{2}(1,idx2(i),:)).*exp(-1j*slope(1)*sub');
                intercept1{ch}(i)=mean(h1.*h2);

            else
                h2=squeeze(channels{2}(2,idx2(i),:)).*exp(-1j*slope(1)*sub');
                intercept1{ch}(i)=mean(h1.*h2);
            end
        end


    end
    
    
    
    mean_intercept=zeros(length(chan),1);   
    for ch=1:1:length(chan)      
       mean_intercept(ch,1)=mean(intercept1{ch});    
    end
    
    % mean_intercept has the channels now

    %% Bartlet: Normal FFT
    figure(1); 
    clf
    subplot(1,3,1);
    plot(freq,unwrap(angle(mean_intercept)),'.');
    H=repmat(mean_intercept.',length(d_vals),1);
    DP=sum(H.*M,2);
    figure(1); 
    subplot(1,3,2);
    plot(d_vals,abs(DP));
    hold on; plot(ground_truth(file_idx)*ones(2,1),[min(abs(DP)),max(abs(DP))],'r');


    %% MUSIC
    h =mean_intercept.';
    H = h'*h;

    [V, D] = eig(H);
    d = d_vals;
    thresh=3;
    nelem = sum(diag(D) > max(diag(D))/thresh);
    for ii=1:length(d)
        e = exp(1i*p_factor*pi*d(ii)./lambda).'; 
        P(ii) = 1./abs(e'*V(:, 1:end-nelem)*V(:, 1:end-nelem)'*e);
    end
    figure(1);
    subplot(1,3,3);
    plot(d, P);
    hold on; plot(ground_truth(file_idx)*ones(2,1),[min(P),max(P)],'r');
    [~,l]=max(P);
    xlabel('Distance');
    ylabel('Power');
    title(num2str(files_all(file_idx)));
    pause(2);
    
    %% Exhaustive
%     options.n_peaks=2;
%     options.res=0.05;
%     options.min_val=0;
%     options.max_val=7.5;
%     options.p_factor=p_factor;
%     
%    
%     if(length(freq)<13)
%        [A,I]=find_sincs_peaks(h.',lambda,options);
%         d1=options.min_val:options.res:options.max_val;
%         h1=figure; hold on; surf(d1,d1,A,'EdgeColor','none');  
%         colormap jet
%     end
    %% Optimization
    options.n_peaks=2;
    options.res=0.05;
    options.min_val=0;
    options.max_val=7.5;
    options.p_factor=p_factor;
    
    npaths=2;
    h1= @(x) optim_fun_lin_distance(x,h.',lambda,options); % Change to optim_fun_lin_distance_approx if you are using 24 frequencies
    problem=createOptimProblem('fmincon','objective',h1,'x0',ones(npaths,1),'lb',zeros(npaths,1),'ub',7.5*ones(npaths,1));
    gs=GlobalSearch;
    run(gs,problem)
    %% Linear optimization
%     options.min_val=0;
%     options.max_val=5;
%     options.alpha=max(abs(h))*0.95;
%     d1=options.min_val:options.res:options.max_val;
%     [v,f]=linear_optim(h,freq,options);
%     figure; hold on; surf(d1,d1,f,'EdgeColor','none');
%     colormap jet
%     figure; hold on; surf(d1,d1,1./f,'EdgeColor','none');
%     colormap jet
    
end

