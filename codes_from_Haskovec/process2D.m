%loads the results of MonteCarlo_2D from 'data2D_K.mat'
%evaluates the statistics and produces the plots

clearvars

%parameter settings for the DBSCAN method
epsilon = 0.05;
minpts = 12;

%max number of memory layers
nK = 8;

%pre-allocation
clsizes = double.empty(nK,0);

for K=1:nK

    %load data
    filename = sprintf('data2D_%d.mat', K);
    load(filename);

    %initialize vdataColl to empty array
    vdataColl = [];
    
    for i=1:nMC
        
        %concatenate vdataColl with vdata
        vdataColl = [vdataColl, vdata];

        %reshape the data into an array of size (2,N)
        xx = reshape(x(i,K,:,:),2,N);

        %calculate distance over the torus
        dist = torusDistances(xx);
        
        %identify clusters using the DBSCAN method with distance 'dist'
        idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

        %plot
        gscatter(xx(1,:),xx(2,:),idx);
        getframe;

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
    
    %process and plot vdataColl
    N = size(vdataColl,1);
    vdataColl(vdataColl==0)=NaN;
    
    mvdata=mean(vdataColl', 'omitnan');
    stvdata=std(vdataColl', 'omitnan');
    
    F = exp( - (1:N) / (pi*int_r*int_r*N) );
    errorbar(F, mvdata, stvdata, 'o','LineWidth',1.5);
    % Create ylabel
    ylabel({'|y^1|'});
    % Create xlabel
    xlabel({'G(\vartheta)'});
    % Create title
    tit = sprintf('K = %d', K);
    title(tit);
    getframe; pause;

    
    %process cluster sizes
    cls = clsizes(K,clsizes(K,:)>0);
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

