function [S, S_low, clusters, Gs] = rate_distortion(G, heuristic, num_pairs)
% Function to compute rate-distortion curve when compressing a given graph
% G.
%
% Input: NxN adjacency matrix G of a possibly weighted, directed network.
% The variable 'heuristic' determines how we choose the pairs of nodes to
% try combining. The variable 'num_pairs' determines how many pairs of
% nodes to try combining at each iteration.
%
% Output: 1xN vector S, where S(n) is the mutual information between a
% random walk on the original network and the same sequence after
% clustering into n clusters. To be precise, S(n) is an upper bound on the
% true mutual information given by assuming that the clustered random walk
% is Markovian. We also return a lower bound S_low on the mutual
% information, the list of clusters, and the network adjacency matrices Gs.
%
% NOTE: Beginning with the full network, the algorithm works by iteratively
% clustering pairs of nodes so as to minimize the upper bound on the mutual
% information. To speed up the algorithm, we only consider pairs of
% clusters specified by 'heuristic'.

% Size of network:
N = size(G,1);
E = sum(G(:))/2;

% Variables to save:
S = zeros(1, N); % Upper bound on entropy rate after clustering
S_low = zeros(1, N); % Lower bound on entropy rate
clusters = cell(1, N); % clusters{n} lists the nodes in each of the n clusters
Gs = cell(1, N); % Gs{n} is the joint transition probability matrix for n clusters

% Transition probability matrix:
P_old = G./repmat(sum(G,2),1,N);

% Compute steady-state distribution:
[p_ss, D] = eigs(P_old'); % Works for all networks
[~, ind] = max(diag(D));
p_ss = p_ss(:,ind)/sum(p_ss(:,ind));
% p_ss = sum(G,2)/sum(G(:)); % Only true for undirected networks
p_ss_old = p_ss;

% Caluclate initial entropy:
logP_old = log2(P_old);
logP_old(isinf(logP_old)) = 0;
S_old = -sum(p_ss_old.*sum(P_old.*logP_old,2));
P_joint = P_old.*repmat(p_ss_old, 1, N);
P_low = P_old;

% Record initial values:
S(end) = S_old;
S_low(end) = S_old;
clusters{end} = num2cell(1:N);
Gs{end} = G;

% Loop over the number of clusterings:
for n = (N-1):-1:2
    
%     tic

    % Different sets of node pairs to try:
    if heuristic == 1
        
        % Try combining all pairs:
        pairs = nchoosek(1:(n+1),2);
        I = pairs(:,1);
        J = pairs(:,2);
        
    elseif heuristic == 2
        
        % Pick num_pairs node pairs at random:
        pairs = nchoosek(1:(n+1),2);
        inds = randsample(size(pairs,1), min([num_pairs, size(pairs,1)])); 
        I = pairs(inds,1);
        J = pairs(inds,2);
        
    elseif heuristic == 3
        
        % Try combining all pairs connected by an edge:
        [I, J] = find(triu(P_old + P_old', 1));
        
    elseif heuristic == 4
        
        % Pick num_pairs node pairs at random that are connected by an edge:
        [I, J] = find(triu(P_old + P_old', 1));
        pair_inds = randsample(length(I), min([num_pairs, length(I)]));
        I = I(pair_inds);
        J = J(pair_inds);
        
    elseif heuristic == 5
        
        % Pick num_pairs node pairs with largest joint transition
        % probabilities:
        P_joint_symm = triu(P_joint + P_joint', 1);
        [~, inds] = maxk(P_joint_symm(:), min([num_pairs, sum(P_joint_symm(:) > 0)]));
        [I, J] = ind2sub([n+1, n+1], inds);
        
    elseif heuristic == 6
        
        % Pick num_pairs node pairs with largest joint transition
        % probabilities plus self-transition probabilities:
        P_joint_symm = triu(P_joint + P_joint' + repmat(diag(P_joint), 1, n+1) +...
            repmat(diag(P_joint)', n+1, 1), 1);
        [~, inds] = maxk(P_joint_symm(:), min([num_pairs, sum(P_joint_symm(:) > 0)]));
        [I, J] = ind2sub([n+1, n+1], inds);
        
    elseif heuristic == 7
        
        % Pick num_pairs node pairs with largest combined stationary
        % probabilities:
        P_ss_temp = triu(repmat(p_ss_old, 1, n+1) + repmat(p_ss_old', n+1, 1), 1);
        [~, inds] = maxk(P_ss_temp(:), min([num_pairs, nchoosek(n+1,2)]));
        [I, J] = ind2sub([n+1, n+1], inds);

    elseif heuristic == 8
        
        % Iteratively add random nodes to one large cluster:
%         I = randsample(n, 1);
        I = 1;
        J = n+1;
        
    else
        error('Variable "setting" is not properly defined.');
    end
    
    % Number of pairs:
    num_pairs_temp = length(I);
    
    % Keep track of all entropies:
    S_all = zeros(1,num_pairs_temp);
    
    % Loop over the pairs of nodes:
    for ind = 1:num_pairs_temp
        
        i = I(ind);
        j = J(ind);
        
        inds_not_ij = [1:(i-1),(i+1):(j-1),(j+1):(n+1)];
        
        % Compute new stationary distribution:
        p_ss_temp = [p_ss_old(inds_not_ij); p_ss_old(i) + p_ss_old(j)];
        
        % Compute new transition probabilities:
        P_temp_1 = sum(repmat(p_ss_old(inds_not_ij), 1, 2).*P_old(inds_not_ij,[i j]), 2);
        P_temp_1 = P_temp_1./p_ss_temp(1:(end-1));
        P_temp_2 = sum(repmat(p_ss_old([i j]), 1, n-1).*P_old([i j], inds_not_ij), 1);
        P_temp_2 = P_temp_2/p_ss_temp(end);
        P_temp_3 = sum(sum(repmat(p_ss_old([i j]), 1, 2).*P_old([i j], [i j])));
        P_temp_3 = P_temp_3/p_ss_temp(end);
        
        logP_temp_1 = log2(P_temp_1);
        logP_temp_1(isinf(logP_temp_1)) = 0;
        logP_temp_2 = log2(P_temp_2);
        logP_temp_2(isinf(logP_temp_2)) = 0;
        logP_temp_3 = log2(P_temp_3);
        logP_temp_3(isinf(logP_temp_3)) = 0;
        
        % Compute change in upper bound on mutual information:
        dS = -sum(p_ss_temp(1:(end-1)).*P_temp_1.*logP_temp_1) - p_ss_temp(end)*sum(P_temp_2.*logP_temp_2) -...
            p_ss_temp(end)*P_temp_3*logP_temp_3 +...
            sum(p_ss_old.*P_old(:,i).*logP_old(:,i)) + sum(p_ss_old.*P_old(:,j).*logP_old(:,j)) +...
            p_ss_old(i)*sum(P_old(i,:).*logP_old(i,:)) + p_ss_old(j)*sum(P_old(j,:).*logP_old(j,:)) -...
            p_ss_old(i)*(P_old(i,i)*logP_old(i,i) + P_old(i,j)*logP_old(i,j)) -...
            p_ss_old(j)*(P_old(j,j)*logP_old(j,j) + P_old(j,i)*logP_old(j,i));
        
        S_temp = S_old + dS;
        
        % Keep track of all entropies:
        S_all(ind) = S_temp;
            
    end
    
    % Find minimum entropy:
    min_inds = find(S_all == min(S_all));
    min_ind = datasample(min_inds, 1);
    
    % Save mutual information:
    S_old = S_all(min_ind);
    S(n) = S_old;
    
    % Compute old transition probabilities:
    i_new = I(min_ind);
    j_new = J(min_ind);
    
    inds_not_ij = [1:(i_new-1),(i_new+1):(j_new-1),(j_new+1):(n+1)];
    
    p_ss_new = [p_ss_old(inds_not_ij); p_ss_old(i_new) + p_ss_old(j_new)];
    
    P_joint = repmat(p_ss_old, 1, n+1).*P_old;
    P_joint = [P_joint(inds_not_ij, inds_not_ij), sum(P_joint(inds_not_ij, [i_new j_new]),2);...
        sum(P_joint([i_new j_new], inds_not_ij),1), sum(sum(P_joint([i_new j_new], [i_new j_new])))];
    P_old = P_joint./repmat(p_ss_new, 1, n);
    p_ss_old = p_ss_new;
    
    logP_old = log2(P_old);
    logP_old(isinf(logP_old)) = 0;
    
    % Record clusters and graph:
    clusters{n} = [clusters{n+1}([1:(i_new-1),(i_new+1):(j_new-1),(j_new+1):(n+1)]),...
        [clusters{n+1}{i_new}, clusters{n+1}{j_new}]];
    Gs{n} = P_joint*2*E;
    
    % Compute lower bound on mutual information:
    P_low = [P_low(:, [1:(i_new-1),(i_new+1):(j_new-1),(j_new+1):(n+1)]),...
        P_low(:,i_new) + P_low(:,j_new)];
    logP_low = log2(P_low);
    logP_low(isinf(logP_low)) = 0;
    S_low(n) = -sum(p_ss.*sum(P_low.*logP_low,2));
    
%     N-n
%     toc
    
end




