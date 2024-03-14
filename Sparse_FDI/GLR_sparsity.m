function [Max_GLR,optimal_GLR_support ] = GLR_sparsity(Z_meas_hist,Rz_Inv,Max_sparsity)%isreal check

[ num_samples,Sample_size] =size(Z_meas_hist);
m = Sample_size;
j1 = 1;
Z_k_n = Z_meas_hist(j1:num_samples,:);
Z_sum = sum(Z_k_n,1);
Max_LLR = -inf;

for loop1 = 1:Max_sparsity
    spasity = loop1;
    %spasity=Max_sparsity;
    Combinations = nchoosek(1:m,spasity);
    [rows_comb,k] = size(Combinations);
    
    Opt_sol = -inf;
    LLr =  0;
    for i1 = 1:rows_comb
        k1 = Combinations(i1,:);
        R = Rz_Inv(k1,k1);
        D = Rz_Inv(:,k1);
        A = (D/R)*D';
        Lamda_kn = 0.5*Z_sum*A*Z_sum'/length(j1:num_samples);
        
        if  ( Lamda_kn > Opt_sol )
            Opt_sol = Lamda_kn;
            LLr = Lamda_kn;
            opt_comb = Combinations(i1,:);
        end
    end
    
    if LLr > Max_LLR
        Max_LLR = LLr;
        detected_support = opt_comb;
    end
    
end
Max_GLR = Max_LLR;
optimal_GLR_support =  detected_support;


