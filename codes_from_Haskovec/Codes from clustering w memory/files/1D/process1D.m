%loads the results of MonteCarlo_1D from 'data1D.mat'
%evaluates the statistics and produces the plots

clearvars

load('data1D.mat');

%parameter settings for the DBSCAN method
epsilon = 0.025;
minpts = 20;

%pre-allocation
clsizes = zeros(nK,1000);
clsizescount = zeros(1,nK);

for K=1:nK
    
    for i=1:nMC

        %reshape the x-data into a vector of length N
        x = reshape(xdata(K,i,:),1,N);

        %calculate distance over the torus
        di = pdist(x');
        di = min(di,1-di);
        dist = squareform(di);

        %calculate particle density
        rho = sum(dist<=int_r) / (int_r*N);
        
        %identify clusters using the DBSCAN method with distance 'dist'
        idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

        %plot
        scatter(x,rho,25,idx,'filled'); axis([0 1 0 max(rho)+1]);
        xlabel('$x$','Interpreter','latex');
        ylabel('$\varrho$','Interpreter','latex');
        getframe;
        %pause;

        %calculate the number of outliers
        outliers(K,i) = sum(idx == -1);

        %calculate the number of clusters
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
            clsizescount(K) = clsizescount(K)+1;
            clsizes(K,clsizescount(K)) = s;
        end

    end

    %process cluster sizes
    cls = clsizes(K,clsizes(K,:)>0);
    if clsizescount(K)==0
        cls=0;
    end
    cls_mean(K) = mean(cls);
    cls_min(K) = double(min(cls));
    cls_max(K) = double(max(cls));

end


%plot cluster sizes
subplot(2,2,1);
errorbar(1:nK, cls_mean, cls_mean-cls_min, -cls_mean+cls_max, 'LineWidth',2);
axis([1 nK -1 400]);
% Create ylabel
ylabel({'cluster size'});
% Create xlabel
xlabel({'K'});
% Create title
title({'min, avg and max cluster size'});

%plot number of clusters
subplot(2,2,2);
errorbar(1:nK, mean(numclusters'), mean(numclusters')-min(numclusters'), -mean(numclusters')+max(numclusters'), 'LineWidth',2);
axis([1 nK -1 16]);
% Create ylabel
ylabel({'number of clusters'});
% Create xlabel
xlabel({'K'});
% Create title
title({'min, avg and max number of clusters'});


%plot number of outliers
subplot(2,2,3);
errorbar(1:nK, mean(outliers'), mean(outliers')-min(outliers'), -mean(outliers')+max(outliers'), 'LineWidth',2);
axis([1 nK -1 400]);
% Create ylabel
ylabel({'number of outliers'});
% Create xlabel
xlabel({'K'});
% Create title
title({'min, avg and max number of outliers'});


%plot percentage of cluster-free outcomes
subplot(2,2,4);
plot(1:nK, 100*sum(numclusters==0,2)/nMC, 'o');
axis([1 nK -10 110]);
% Create ylabel
ylabel({'percentage'});
% Create xlabel
xlabel({'K'});
% Create title
title({'% of cluster-free outcomes'});
