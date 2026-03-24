%loads the results of MonteCarlo_2D from 'data2D_K.mat'
%evaluates the statistics and produces the plots

clearvars

%parameter settings for the DBSCAN method
epsilon = 0.05;
minpts = 12;

%Number of runs in each Monte-Carlo simulation
nMC = 100;

%number of values of tau
taus = 0:30:300;
Ntaus = size(taus,2);
delayType = "Reaction";
d = 2;
fileName = "MCData";

%number of particles
N=400;

%pre-allocation
clsizes = double.empty(Ntaus,0);

for K=1:Ntaus

    tau = taus(K);

    %load data
    folderPath = MCFolderPath(delayType,d,tau);
    postfix = MCFilePostfix(delayType,d,tau);
    filePath = MCFilePath(folderPath,fileName,postfix);
    load(filePath);
    
    for i=1:nMC

        %reshape the data into an array of size (2,N)
        xx = results{i}.xRec(:,:,end).';

        %calculate distance over the torus
        dist = torusDistances(xx.');
        
        %identify clusters using the DBSCAN method with distance 'dist'
        idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

        %plot
        gscatter(xx(1,:),xx(2,:),idx); getframe;
        %pause;

        %number of outliers
        outliers(K,i) = sum(idx == -1);

        %number of clusters
        numclusters(K,i) = size(unique(idx(idx>0)),1);

        %stdev of cluster sizes
        Midx=idx*ones(1,numclusters(K,i));
        Mclu=(1:numclusters(K,i))'*ones(1,N);
        clustsizes = sum(Midx==Mclu');
        clstss(K,i,1:N) = [clustsizes, -ones(1,N-numel(clustsizes))];
        stdev_clusizes(K,i) = std(clustsizes);

        %cluster sizes
        for j=1:numclusters(K,i)
            s = sum(idx==j);
            clsizes(K,end+1) = s;
        end

    end

    
    %process cluster sizes
    cls = clsizes(K,clsizes(K,:)>0);
    cls_mean(K) = mean(cls);
    cls_min(K) = double(min(cls));
    cls_max(K) = double(max(cls));

end



%plot cluster sizes
subplot(2,2,1);
errorbar(1:Ntaus, cls_mean, cls_mean-cls_min, -cls_mean+cls_max, 'LineWidth',2);
axis([1 Ntaus -1 400]);
% Create ylabel
ylabel({'cluster size'});
% Create xlabel
xlabel({'tau'});
% Create title
title({'min, avg and max cluster size'});

%plot number of clusters
subplot(2,2,2);
errorbar(1:Ntaus, mean(numclusters'), mean(numclusters')-min(numclusters'), -mean(numclusters')+max(numclusters'), 'LineWidth',2);
axis([1 Ntaus -1 16]);
% Create ylabel
ylabel({'number of clusters'});
% Create xlabel
xlabel({'K'});
% Create title
title({'min, avg and max number of clusters'});


%plot number of outliers
subplot(2,2,3);
errorbar(1:Ntaus, mean(outliers'), mean(outliers')-min(outliers'), -mean(outliers')+max(outliers'), 'LineWidth',2);
axis([1 Ntaus -1 400]);
% Create ylabel
ylabel({'number of outliers'});
% Create xlabel
xlabel({'K'});
% Create title
title({'min, avg and max number of outliers'});


%plot percentage of cluster-free outcomes
subplot(2,2,4);
plot(1:Ntaus, 100*sum(numclusters==0,2)/nMC, 'o');
axis([1 Ntaus -10 110]);
% Create ylabel
ylabel({'percentage'});
% Create xlabel
xlabel({'K'});
% Create title
title({'% of cluster-free outcomes'});

