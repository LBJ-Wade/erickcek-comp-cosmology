%this one finds top100 halos and their subhalos for each
halos = csvread('T30MC40.ascii');
sorted_halos = sortrows(halos,3); %sorts by third column ie mass!!
sorted_halos = flipud(sorted_halos);
top100 = sorted_halos(1:100,:);
% sub_halos = np.zeros(100);
% for i=1:length(sorted_halos)
%     host = halos(i,13);
%     if (host == -1)
%     else
%         sub_halos(host) = sub_halos(host) + 1;
%     end
%         
% end