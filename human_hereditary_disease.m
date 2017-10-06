% Loading the dataset 
Data = importdata('/rmt/csfiles/pgrads/mbva620/3d_filtered_outfile', ' ');
% Loading the disease types 
labels = importdata('/rmt/csfiles/pgrads/mbva620/labels', ' ');
label = zeros(length(labels),1);
% Filter the unique set of labels
unique_labels = unique(labels);
% Change the datatype of label 
for i = 1:length(unique_labels)
   for j = 1:length(labels)
       if isequal(labels(j,1),unique_labels(i,1))
           label(j,1) = i;
       end
   end
end
[r,c] = size(Data);
% Number of clusters required K value
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
% Cluster groups are initialised   
    cluster1 = zeros(K,length(Data));
    cluster2 = zeros(K,length(Data));
% Iteration count is intialised
    Iteration_count = 0;
    % Flag variable
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
        
        %find the index of the cluster where distance is minimum
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
                % Save the clusters in a text file 'disease_cluster.txt'
                save disease_cluster.txt Cluster_Matrix -ascii -append;
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
        
     figure(1);
    
        scatter3(Data(:,1),Data(:,2),Data(:,3),50,label(:,1),'filled');
        title('Disease plot - Each color represents a disease type');
      
    figure(2);
    Color_map=colormap;
    vec= 1:floor(64/K):64;
    for i = 1:K
    scatter3(Data(x == i,1),Data(x == i,2),Data(x == i,3),50,Color_map(vec(i),:),'filled');
    title('Cluster groups - Each color represents a cluster group');
    hold on;
    end
    SSD = 0;
    % Find the squared distance
    if choice == 1
        SD = distance;
    elseif choice == 2
        SD = distance .^ 2;
    end
    
    for i = 1:K
        % Find the sum squared distance
        SSD = SSD + sum(SD(i,(x == i)));
    end
        % Sum squared distance vector
        SumSquaredDist(Iteration_count) = SSD;
    end

figure(3)
plot(SumSquaredDist); 
title(' Sum of Squared distance Vs No of Iterations');
xlabel('No.of Iterations');
ylabel('Sum Squared Distance');
