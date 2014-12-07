function [channels, timestamp,packet_ids,MACs,sec,usec]=process_trace_channel(filename,plot_channels)

intel_log = read_bf_file_channel(filename);
num_entries = length(intel_log);
disp(num_entries);
csi = cell(1,num_entries);
MACs=zeros(1,num_entries);
packet_ids=zeros(1,num_entries);
channels=zeros(2,num_entries,30);
num_tx = zeros(1,num_entries);
timestamp = zeros(1,num_entries);
sec=zeros(1,num_entries);
usec=zeros(1,num_entries);
ii = 0;
for i = 1:num_entries
    if(~isempty(intel_log{i}))
      %  if (isempty(find(intel_log{i}.csi(1,1:2,:)==0, 1)))
            ii = ii+1;
            timestamp(ii) = intel_log{i}.timestamp_low;
            MACs(ii) = intel_log{i}.mac;
            packet_ids(ii) = intel_log{i}.packet_id(end);
            csi{ii} = intel_log{i}.csi;
            channels(1:2,ii,:) = reshape(csi{ii}(1,1:2,:),2,30);
            num_tx(ii) = intel_log{i}.Ntx;
            sec(ii)=intel_log{i}.tv_sec;
            usec(ii)=intel_log{i}.tv_usec;
       % else
        %    fprintf(1,'here');
       % end
    else
        fprintf(1,'here');
    end
end
cutoff = find(timestamp>0,1,'last');
timestamp = timestamp(1:cutoff);
MACs = MACs(1:cutoff);
packet_ids = packet_ids(1:cutoff);
channels = channels(:,1:cutoff,:);
sec=sec(1:cutoff);
usec=usec(1:cutoff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (intel_log{1}.perm(1)==2)
%     %channels([1,2],:,ch)=channels([2,1],:,ch);
%     disp('switch rx1 and rx2');
% end
% if (intel_log{1}.perm(1)==3||intel_log{1}.perm(2)==3)
%     fprintf(1,'\n*****ANTENNA 3 POTENTIALLY HAS HIGH SNR*****\n');
% end

phones = unique(MACs);
for ch = plot_channels
    figure;
%    subplot(1,2,1);
    plot(timestamp(MACs==phones(1))/1e6,unwrap(angle(channels(1,(MACs==phones(1)),ch)./channels(2,(MACs==phones(1)),ch))));
%    subplot(1,2,2);
%    plot(timestamp(MACs==phones(2))/1e6,unwrap(angle(channels(1,(MACs==phones(2)),ch)./channels(2,(MACs==phones(2)),ch))));
    
end
%figure; plot(timestamp/1e6,ones(size(timestamp)),'.');
