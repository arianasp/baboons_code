function [ state_mat ] = generate_state_mat( N )
%Generates a matrix of all possible states, for easy multiplication later on.
%In this list, each row represents one state, and the numbers given are
%what is going to be multipled by the alpha and beta for each state. The
%rows look like [x11 x12 x13 ... x21 x22 x23 ... xN1 xN2 xN3 ... xNN];
%Note that xii = NaN
%INPUTS:
%   N: number of individuals / spins
%OUTPUTS:
%   state_mat: [(N + N^2) x 2^N matrix] as described above

state_mat = zeros(N+(N^2-N)/2,2^N);

for idx = 1:(2^N)
    xi = dec2bin(idx-1,N) - '0'; %convert decimal index to a binary vector 
    
    xi(find(xi==0)) = -1; %convert 0's to -1's for computation
    
    xij = xi' * xi; %generate couplings matrix
    
    tri_mat = triu(ones(N,N),1); %generate upper triangular matrix to use for indexing into couplings matrix
    
    xij = xij(find(tri_mat==1)); %get upper triangular part of xij matrix as a vector
    
    state_mat(:,idx) = [xi xij']'; %add column to state matrix
    
end


end

