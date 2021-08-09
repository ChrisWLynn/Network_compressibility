% Script to compute rate-distortion curves and other measurements for
% directed samples of real networks
%
% NOTE: The difference between this script and
% 'script_rate_distortion_real_samples' is that here we consider directed
% versions of the networks
%
% NOTE: We use this for the analysis of directed networks in the paper.

% Graph info:
info_struct = load('graphs_info.mat');
info = info_struct.graph_info;

% Limit to directed networks:
info = info(find(info.directed),:);

% Graph families:
families = unique(info.family);
num_families = length(families);

% Graphs to use top joint transition probability heuristic:
graphs_topTrans = {'G_FOLDOC', 'G_Clavier', 'G_SC_2'};

% Loop over graph families:
for i = 1:num_families
    
    % Load graphs:
    graph_struct = load(['graphs_', families{i}, '_samples_directed.mat']);
    
    % Different graphs in this family:
    graph_names = fieldnames(graph_struct);
    num_graphs = length(graph_names);
    
    % Different things to compute:
    rate_distortion_upper = cell(num_graphs, 1);
    rate_distortion_lower = cell(num_graphs, 1);
    avg_outDeg_inside = cell(num_graphs, 1);
    avg_inDeg_inside = cell(num_graphs, 1);
    avg_outDeg_outside = cell(num_graphs, 1);
    avg_inDeg_outside = cell(num_graphs, 1);
    num_edges_inside = cell(num_graphs, 1);
    num_edges_boundary_out = cell(num_graphs, 1);
    num_edges_boundary_in = cell(num_graphs, 1);
    num_edges_outside = cell(num_graphs, 1);
    cluster_size = cell(num_graphs, 1);
    compressibility = cell(num_graphs, 1);
    avg_deg = cell(num_graphs, 1);
    avg_clustering = cell(num_graphs, 1);
    outDeg_het = cell(num_graphs, 1);
    inDeg_het = cell(num_graphs, 1);
    
    % Structure to save:
    measures_struct = [];
    
    % Loop over different graphs in this family:
    for j = 1:num_graphs
        
        % Load graph samples:
        Gs = graph_struct.(graph_names{j});
        num_samples = length(Gs);
        
        % Things to compute for this graph:
        rate_distortion_upper_temp = cell(num_samples,1);
        rate_distortion_lower_temp = cell(num_samples,1);
        avg_outDeg_inside_temp = cell(num_samples,1);
        avg_inDeg_inside_temp = cell(num_samples,1);
        avg_outDeg_outside_temp = cell(num_samples,1);
        avg_inDeg_outside_temp = cell(num_samples,1);
        num_edges_inside_temp = cell(num_samples,1);
        num_edges_boundary_out_temp = cell(num_samples,1);
        num_edges_boundary_in_temp = cell(num_samples,1);
        num_edges_outside_temp = cell(num_samples,1);
        cluster_size_temp = cell(num_samples,1);
        compressibility_temp = zeros(num_samples, 1);
        avg_deg_temp = zeros(num_samples, 1);
        avg_clustering_temp = zeros(num_samples, 1);
        outDeg_het_temp = zeros(num_samples, 1);
        inDeg_het_temp = zeros(num_samples, 1);
        
        % Loop over graph samples:
        for k = 1:num_samples
            
            tic
            
            % Load graph sample:
            G = Gs{k};
            N = size(G,1);
            ks_out = sum(G,2);
            ks_in = sum(G,1);
        
            % Compute simple graph measures:
            avg_deg_temp(k) = mean(ks_out);
            avg_clustering_temp(k) = avgClusteringCoefficient(double((G+G')>0));
            outDeg_het_temp(k) = sum(sum(abs(ks_out - ks_out')))/(N*(N-1)*mean(ks_out));
            inDeg_het_temp(k) = sum(sum(abs(ks_in - ks_in')))/(N*(N-1)*mean(ks_in));
            
            % Compute rate-distortion curve:
            if ismember(graph_names{j}, graphs_topTrans)
                
                % Use top joint transition probability pairs heuristic:
                [rd_upper, rd_lower, Cs, ~] = rate_distortion(G, 5, 100);
                
            else
                
                % Use top stationary distribution pairs heuristic:
                [rd_upper, rd_lower, Cs, ~] = rate_distortion(G, 7, 100);
                
            end
        
            rate_distortion_upper_temp{k} = rd_upper;
            rate_distortion_lower_temp{k} = rd_lower;
            compressibility_temp(k) = mean(rd_upper(end) - rd_upper);
        
            % Compute statistics about large cluster:
            avg_outDeg_inside_temp_temp = zeros(1, N-2);
            avg_inDeg_inside_temp_temp = zeros(1, N-2);
            avg_outDeg_outside_temp_temp = zeros(1, N-2);
            avg_inDeg_outside_temp_temp = zeros(1, N-2);
            num_edges_inside_temp_temp = zeros(1, N-1);
            num_edges_boundary_out_temp_temp = zeros(1, N-1);
            num_edges_boundary_in_temp_temp = zeros(1, N-1);
            num_edges_outside_temp_temp = zeros(1, N-1);
            cluster_size_temp_temp = zeros(1, N);
            
            cluster_size_temp_temp(1) = N;
            cluster_size_temp_temp(end) = 1;
            
            num_edges_inside_temp_temp(1) = sum(sum(G))/2;
            num_edges_boundary_out_temp_temp(1) = 0;
            num_edges_boundary_in_temp_temp(1) = 0;
            num_edges_outside_temp_temp(1) = 0;
        
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
                cluster_size_temp_temp(l) = length(cluster);
                avg_outDeg_inside_temp_temp(l-1) = mean(ks_out(cluster));
                avg_inDeg_inside_temp_temp(l-1) = mean(ks_in(cluster));
                avg_outDeg_outside_temp_temp(l-1) = mean(ks_out(not_cluster));
                avg_inDeg_outside_temp_temp(l-1) = mean(ks_in(not_cluster));
                num_edges_inside_temp_temp(l) = sum(sum(G(cluster, cluster)));
                num_edges_boundary_out_temp_temp(l) = sum(sum(G(cluster, not_cluster)));
                num_edges_boundary_in_temp_temp(l) = sum(sum(G(not_cluster, cluster)));
                num_edges_outside_temp_temp(l) = sum(sum(G(not_cluster, not_cluster)));
                
            end
            
            avg_outDeg_inside_temp{k} = avg_outDeg_inside_temp_temp;
            avg_inDeg_inside_temp{k} = avg_inDeg_inside_temp_temp;
            avg_outDeg_outside_temp{k} = avg_outDeg_outside_temp_temp;
            avg_inDeg_outside_temp{k} = avg_inDeg_outside_temp_temp;
            num_edges_inside_temp{k} = num_edges_inside_temp_temp;
            num_edges_boundary_out_temp{k} = num_edges_boundary_out_temp_temp;
            num_edges_boundary_in_temp{k} = num_edges_boundary_in_temp_temp;
            num_edges_outside_temp{k} = num_edges_outside_temp_temp;
            cluster_size_temp{k} = cluster_size_temp_temp;

            i
            j
            k
            toc

        end
        
        % Record measures for this sample:
        rate_distortion_upper{j} =  rate_distortion_upper_temp;
        rate_distortion_lower{j} = rate_distortion_lower_temp;
        avg_outDeg_inside{j} = avg_outDeg_inside_temp;
        avg_inDeg_inside{j} = avg_inDeg_inside_temp;
        avg_outDeg_outside{j} = avg_outDeg_outside_temp;
        avg_inDeg_outside{j} = avg_inDeg_outside_temp;
        num_edges_inside{j} = num_edges_inside_temp;
        num_edges_boundary_out{j} = num_edges_boundary_out_temp;
        num_edges_boundary_in{j} = num_edges_boundary_in_temp;
        num_edges_outside{j} = num_edges_outside_temp;
        cluster_size{j} = cluster_size_temp;
        compressibility{j} = compressibility_temp;
        avg_deg{j} = avg_deg_temp;
        avg_clustering{j} = avg_clustering_temp;
        outDeg_het{j} = outDeg_het_temp;
        inDeg_het{j} = inDeg_het_temp;
        
    end
    
    % Add measures to structure:
    measures_struct.([families{i}, '_rate_distortion_upper']) = rate_distortion_upper;
    measures_struct.([families{i}, '_rate_distortion_lower']) = rate_distortion_lower;
    measures_struct.([families{i}, '_avg_outDeg_inside']) = avg_outDeg_inside;
    measures_struct.([families{i}, '_avg_inDeg_inside']) = avg_inDeg_inside;
    measures_struct.([families{i}, '_avg_outDeg_outside']) = avg_outDeg_outside;
    measures_struct.([families{i}, '_avg_inDeg_outside']) = avg_inDeg_outside;
    measures_struct.([families{i}, '_num_edges_inside']) = num_edges_inside;
    measures_struct.([families{i}, '_num_edges_boundary_out']) = num_edges_boundary_out;
    measures_struct.([families{i}, '_num_edges_boundary_in']) = num_edges_boundary_in;
    measures_struct.([families{i}, '_num_edges_outside']) = num_edges_outside;
    measures_struct.([families{i}, '_cluster_size']) = cluster_size;
    measures_struct.([families{i}, '_compressibility']) = compressibility;
    measures_struct.([families{i}, '_avg_deg']) = avg_deg;
    measures_struct.([families{i}, '_avg_clustering']) = avg_clustering;
    measures_struct.([families{i}, '_outDeg_het']) = outDeg_het;
    measures_struct.([families{i}, '_inDeg_het']) = inDeg_het;
    
    % Save measures:
    save(['graphs_', families{i}, '_rate_distortion_samples_directed'], '-struct', 'measures_struct');
    
end


