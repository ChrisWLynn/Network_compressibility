% Script to compute rate-distortion curves and other measurements for
% undirected samples of real networks
%
% NOTE: We use this for the main analysis in the paper.

% Graph info:
info_struct = load('graphs_info.mat');
info = info_struct.graph_info;

% Graph families:
families = unique(info.family);
num_families = length(families);

% Graphs to use top joint transition probability heuristic:
graphs_topTrans = {'G_FOLDOC', 'G_AstroPh', 'G_Facebook', 'G_karate',...
    'G_Clavier', 'G_SC_2', 'G_macaques', 'G_zebras'};

% Loop over graph families:
for i = 1:num_families
    
    % Load graphs:
    graph_struct = load(['graphs_', families{i}, '_samples_undirected.mat']);
    
    % Different graphs in this family:
    graph_names = fieldnames(graph_struct);
    num_graphs = length(graph_names);
    
    % Different things to compute:
    rate_distortion_upper = cell(num_graphs, 1);
    rate_distortion_lower = cell(num_graphs, 1);
    avg_deg_inside = cell(num_graphs, 1);
    avg_deg_outside = cell(num_graphs, 1);
    num_edges_inside = cell(num_graphs, 1);
    num_edges_boundary = cell(num_graphs, 1);
    num_edges_outside = cell(num_graphs, 1);
    cluster_size = cell(num_graphs, 1);
    compressibility = cell(num_graphs, 1);
    deg_het = cell(num_graphs, 1);
    avg_clustering = cell(num_graphs, 1);
    avg_deg = cell(num_graphs, 1);
    
    % Structure to save:
    measures_struct = [];
    
    % Loop over different graphs in this family:
    for j = 1:num_graphs
        
        % Load graph samples:
        Gs = graph_struct.(graph_names{j});
        num_samples = length(Gs);
        N = size(Gs{1},1);
        
        % Things to compute for this graph:
        rate_distortion_upper_temp = zeros(num_samples, N);
        rate_distortion_lower_temp = zeros(num_samples, N);
        avg_deg_inside_temp = zeros(num_samples, N-2);
        avg_deg_outside_temp = zeros(num_samples, N-2);
        num_edges_inside_temp = zeros(num_samples, N-1);
        num_edges_boundary_temp = zeros(num_samples, N-1);
        num_edges_outside_temp = zeros(num_samples, N-1);
        cluster_size_temp = zeros(num_samples, N);
        compressibility_temp = zeros(num_samples, 1);
        deg_het_temp = zeros(num_samples, 1);
        avg_clustering_temp = zeros(num_samples, 1);
        avg_deg_temp = zeros(num_samples, 1);
        
        % Loop over graph samples:
        for k = 1:num_samples
            
            tic
            
            % Load graph sample:
            G = Gs{k};
            ks = sum(G,1);
        
            % Compute simple graph measures:
            deg_het_temp(k) = sum(sum(abs(ks - ks')))/(N*(N-1)*mean(ks));
            avg_clustering_temp(k) = avgClusteringCoefficient(G);
            avg_deg_temp(k) = mean(ks);
            
            % Compute rate-distortion curve:
            if ismember(graph_names{j}, graphs_topTrans)
                
                % Use top joint transition probability pairs heuristic:
                [rd_upper, rd_lower, Cs, ~] = rate_distortion_upper(G, 5, 100);
                
            else
                
                % Use top stationary distribution pairs heuristic:
                [rd_upper, rd_lower, Cs, ~] = rate_distortion_upper(G, 7, 100);
                
            end
        
            rate_distortion_upper_temp(k,:) = rd_upper;
            rate_distortion_lower_temp(k,:) = rd_lower;
            compressibility_temp(k) = mean(rd_upper(end) - rd_upper);
        
            % Compute statistics about large cluster:
            cluster_size_temp(k,1) = N;
            cluster_size_temp(k,end) = 1;
            
            num_edges_inside_temp(k,1) = sum(sum(G))/2;
            num_edges_boundary_temp(k,1) = 0;
            num_edges_outside_temp(k,1) = 0;
        
            for l = 2:(N-1)
                
                C = Cs{l};
                
                % Get largest cluster:
                sizes = zeros(1, length(C));
                for m = 1:length(C)
                    sizes(m) = length(C{m});
                end
                [~, ind] = max(sizes);
                cluster = C{ind};
                not_cluster = setdiff(1:N, cluster);
                
                % Compute things:
                cluster_size_temp(k,l) = length(cluster);
                avg_deg_inside_temp(k,l-1) = mean(ks(cluster));
                avg_deg_outside_temp(k,l-1) = mean(ks(not_cluster));
                num_edges_inside_temp(k,l) = sum(sum(G(cluster, cluster)))/2;
                num_edges_boundary_temp(k,l) = sum(sum(G(cluster, not_cluster)));
                num_edges_outside_temp(k,l) = sum(sum(G(not_cluster, not_cluster)))/2;
                
            end

            i
            j
            k
            toc

        end
        
        % Record measures for this sample:
        rate_distortion_upper{j} =  rate_distortion_upper_temp;
        rate_distortion_lower{j} = rate_distortion_lower_temp;
        avg_deg_inside{j} = avg_deg_inside_temp;
        avg_deg_outside{j} = avg_deg_outside_temp;
        num_edges_inside{j} = num_edges_inside_temp;
        num_edges_boundary{j} = num_edges_boundary_temp;
        num_edges_outside{j} = num_edges_outside_temp;
        cluster_size{j} = cluster_size_temp;
        compressibility{j} = compressibility_temp;
        deg_het{j} = deg_het_temp;
        avg_clustering{j} = avg_clustering_temp;
        avg_deg{j} = avg_deg_temp;
        
    end
    
    % Add measures to structure:
    measures_struct.([families{i}, '_rate_distortion_upper']) = rate_distortion_upper;
    measures_struct.([families{i}, '_rate_distortion_lower']) = rate_distortion_lower;
    measures_struct.([families{i}, '_avg_deg_inside']) = avg_deg_inside;
    measures_struct.([families{i}, '_avg_deg_outside']) = avg_deg_outside;
    measures_struct.([families{i}, '_num_edges_inside']) = num_edges_inside;
    measures_struct.([families{i}, '_num_edges_boundary']) = num_edges_boundary;
    measures_struct.([families{i}, '_num_edges_outside']) = num_edges_outside;
    measures_struct.([families{i}, '_cluster_size']) = cluster_size;
    measures_struct.([families{i}, '_compressibility']) = compressibility;
    measures_struct.([families{i}, '_deg_het']) = deg_het;
    measures_struct.([families{i}, '_avg_clustering']) = avg_clustering;
    measures_struct.([families{i}, '_avg_deg']) = avg_deg;
    
    % Save measures:
    save(['graphs_', families{i}, '_rate_distortion_samples'], '-struct', 'measures_struct');
    
end


