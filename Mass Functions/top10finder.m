%this script looks at our microhalo populations sorted to find top 10
%largest in each of the 4 512 .000030 Mpc sims ran 
cats = dir('*.ascii');
for i=1:length(cats)
    halos = csvread(cats(i).name);
    sorted_halos = sortrows(halos,3); %sorts by third column ie mass!!
    sorted_halos = flipud(sorted_halos);
    top10 = sorted_halos(1:10,:);
    top10(:,3);
end