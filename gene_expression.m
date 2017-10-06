% Loading the dataset 
data = importdata('/rmt/csfiles/pgrads/mbva620/EcoliDatasetCW5810_truncated.txt', '\t');
% The names of the genes
gene_names = data.textdata;
% Pre-process the data
% Discard column-2 and the first column is gene name
gene_name = gene_names(:,1);
% The numerical gene values
gene_values = data.data;
% Average 18 column values in groups of 3
[rows,cols] = size(gene_values);
r=reshape(gene_values,rows,3,6);
Data=squeeze(mean(r,2)); % The pre-processed data is stored in 'Data'
[r,c] = size(Data);
K = input('Enter the number of clusters: ');
% Choose the distance metric among the 2 different distance metric to cluster the genes
% Choice-1 : Euclidean Distance
% Choice-2 : Correlation   
 choice = input('Enter your choice: 1. Euclidean distance 2. Correlation');
% Choose K random centroids
rand_index = randperm(length(Data),K);
Centroids = Data(rand_index,:);
% Centroid matrices are initialised for forming the stopping criterion
Centroid1 =Centroids;
Centroid2 = zeros(r,c);
% Cluster groups are initialised which helps in forming the cluster matrix
cluster1 = zeros(K,length(Data));
cluster2 = zeros(K,length(Data));
% Initialising the iteration count
    Iteration_count = 0;
% Initialising flag variable
flag = 0;
    while flag ~= 1
      Iteration_count = Iteration_count + 1;  
      % Use distance metrics to calculate the distance 
        if choice == 1
            distance = pdist2(Centroids,Data,'euclidean');
        elseif choice == 2
            distance = pdist2(Centroids,Data,'correlation');
        else
            disp ('Invalid Choice'); 
        end
        
        % find the index of the cluster where distance is minimum
        % The index denotes the cluster in which the data point belongs 
        [min_distance,index] = min(distance);
                  
        % The dimensions of the distance matrix is calculated
        [row,col] = size(distance);

        for i = 1:col
         cluster2(index(i),i) = 1;
        end
            
            for i = 1:K
                % Find the cluster index to form the cluster matrix
                Cluster_Index = find(cluster2(i,:)==1);
                % Cluster matrix
                Cluster_Matrix = zeros(length(Cluster_Index),c);
                for j = 1:length(Cluster_Index)
                    Cluster_Matrix(j,:) = Data(Cluster_Index(j),:);
                end
                % Save the clusters in a text file 'gene_cluster.txt'
                save gene_cluster.txt Cluster_Matrix -ascii -append;
                %update centroid co-ordinates
                Centroids(i,:) = sum(Cluster_Matrix)/length(Cluster_Index);
                Centroid2 = Centroids;
            end
            % Compare the centroids to form a stopping criterion
       if Centroid1 == Centroid2
           flag = 1;
       end
            Centroid1 = Centroid2;
            cluster1 = cluster2;
            cluster2 = zeros(K,length(Data));
            
        v = 1:1:K;
        x = v*cluster1;
        
        for i = 1:K
            figure(i);
            trans_matrix = Data(x==i,:).';
            % Plot the cluster 
            plot(trans_matrix);
            title('Clusters');
        end
    SSD = 0;
    % Find the squared distance
    if choice == 1
        SD = distance;
    elseif choice == 2
        SD = power(distance,2);
    end
    % Find the sum squared distance
    for i = 1:K
        SSD = SSD + sum(SD(i,(x == i)));
    end
        % Sum squared distance vector
        SumSquaredDist(Iteration_count) = SSD;
    end
figure(K+1)
plot(SumSquaredDist);
xlabel('No.of Iterations');
ylabel('Sum Squared Distance');
title(' Sum of Squared distance Vs No of Iterations');