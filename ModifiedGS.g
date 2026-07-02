AreColumnsOrthogonal := function(mat)
    local col_count, i, j, col_i, col_j, dot_product, entry_index, non_orthog_cols_found;
    
    col_count := Length(mat[1]);  # Number of columns

    for i in [1..col_count] do
        col_i := mat{[1..col_count]}{[i]};
        
        
        for j in [i..col_count] do
            if j = i then
                continue;
            fi;
            col_j := mat{[1..col_count]}{[j]};
            
            #Compute the dot product of col_i and col_j
            dot_product := 0;
            for entry_index in [1..Length(col_i)] do
                #Print("Col i entry: ", col_i[entry_index], " Col j entry: ", col_j[entry_index], "\n");
                dot_product := dot_product + ComplexConjugate(col_i[entry_index])*(col_j[entry_index]); #this is the standard complex dot product
            od;
            #-------------------------ALL FOR DEBUGGING--------------
            # Check if dot_product is 0 or close to zero
            # if not IsZero(dot_product) and dot_product*ComplexConjugate(dot_product) > Rat(0.01) then
            #     Print("Cols ", i, " and ", j, " are not orthogonal! \n");
            #     Print("The dot product was: ", dot_product, "\n");
            #     Print("Cols ", i, " and ", j, " are not orthogonal! \n");
            #     Print("The dot product was: ", dot_product, "\n");
            #     non_orthog_cols_found := true;    
            # fi;

            if not IsZero(dot_product) then
                Print("Cols ", i, " and ", j, " are not orthogonal! \n");
                Print("The dot product was: ", dot_product, "\n");
                non_orthog_cols_found := true;
                
            fi;

        od;
    od;
    if non_orthog_cols_found = true then
        return false;
    fi;
    return true;
    
end;

#A function to compute the complex norm squared of a vector
ComplexNormSquared:= function(col)
    local col_norm_squared, entry, entry_norm;
    col_norm_squared := 0;
    for entry in col do
        entry_norm := entry*ComplexConjugate(entry);
        col_norm_squared := col_norm_squared + entry_norm;
    od;
    return col_norm_squared;
end;

#MODIFIED GRAM SCHMIDT ALGORITHM
#Takes as input
# 1) a full rank intertwiner between the tensor product of a pair of irreps (a,b) and their block diagonal form
# 2) a list of distinct fusion outcomes of a & b as strings (with multiplicities).
# Example for point (2)- if a*b = 3c + 2d + e; then the decomposition_list paramter would be ["e", "2d", "3c"]
Modified_Gram_Schmidt := function(mat, decomposition_list)
    #Local Variables
    local
            col_norm_squared, col, col_norm, entry, entry_norm,
            distinct_fusion_outcome_count, intertwiner_dimension,
            orthogonalized_intertwiner, orthogonalized_intertwiner_index,
            col_index, i, n_i, d_c_i, fusion_outcome_label,
            orth_projector_list, multiplicity_index, u_col_index,
            Pi_nu, t_nu_transpose, row, f_mu_r, outer_product,
            f_mu_1, alpha_mu, nu, t_mu, f_mu, beta_mu,
            e_mu, e_mu_list, row_number, col_number, F_i, mu, T_i, P, Pi_mu, f_mu_transpose, T_i_index; 

    #INITIALIZE VARIABLES
    #Record the total number of distinct fusion outcomes
    distinct_fusion_outcome_count := Length(decomposition_list);
    #Record dimension of the input matrix
    intertwiner_dimension := Length(mat);
    #Once the cols of the input matrix "mat" are orthogonalized, the result will be stored in orthogonalized_intertwiner
    orthogonalized_intertwiner := NullMat(intertwiner_dimension, intertwiner_dimension);
    orthogonalized_intertwiner_index := 0; #this index will help populate the orthogonalized_intertwiner matrix later

    col_index := 1; #a way to index the cols of the input matrix "mat" (which is the full rank intertwiner we wish to orthogonalize)

    #MAP T_i -> F_i for i in {1...distinct_fusion_outcome_count}
    for i in [1..distinct_fusion_outcome_count] do
        
        #Initialize necessary information about this fusion outcome
        n_i := label_info(decomposition_list[i])[1]; #the multiplicity of the irrep c_i in the decomposition
        fusion_outcome_label := label_info(decomposition_list[i])[2]; #the label of the irrep c_i in the decomposition
        d_c_i := DegreeOfRepresentation(irrep_dictionary.(fusion_outcome_label)); #the degree of the irrep c_i
        
        #Pick out T_i 
        T_i := mat{[1..intertwiner_dimension]}{[col_index..(col_index + d_c_i*n_i - 1)]}; #Picks out columns of the input matrix between indices col_index to (col_index + d_c_i*n_i) inclusive
        T_i_index := 0; #will be used to index the columns of T_i

        if n_i = 1 then
            #All columns are already mutually orthogonal, set F_i = T_i
            F_i := T_i;
            col_index := col_index + d_c_i;
        else
            #Implement GS1
            #Initialize a list of Orthogonal Projectors Pi_nu for nu = {1..n_i - 1}
            P:= IdentityMat(intertwiner_dimension);
            #T_i is mapped to F_i via orthogonalization
            F_i := NullMat(intertwiner_dimension, n_i*d_c_i); 
            for mu in [1..n_i] do
                #Define t_mu (t_mu is the mu^{th} block of d_c_i cols of T_i)
                t_mu := mat{[1..intertwiner_dimension]}{[col_index..(col_index + d_c_i - 1)]}; #this extracts columns from col_index to (col_index+d_c_i) with endpoints inclusive;
                
                #f_mu is the orthogonalized version of t_mu
                f_mu := t_mu; #initialize f_mu as t_mu
                
                if mu = 1 then
                    f_mu := t_mu; #let f_mu stay as t_mu
                else
                    f_mu := P*t_mu;
                fi;

                

                #---------------Compute Pi_mu and add it to the orth_projector_list---------------
                #Update P
                Pi_mu := NullMat(intertwiner_dimension, intertwiner_dimension);
                
                #(1) Compute the projector onto the col span of {f_mu}
                #We must iterate over cols of f_mu to calculate Pi_nu, but indexing rows is easier
                f_mu_transpose := TransposedMat(f_mu); #indexing rows is much easier than cols
                for row in f_mu_transpose do
                    f_mu_r := TransposedMat([row]); #This is the rth column of f_mu (r ranges from 1 to d_c_i in this loop)
                    outer_product := f_mu_r*ComplexConjugate([row]) ;#The orthogonal projector onto t_nu_r
                    Pi_mu := Pi_mu + outer_product ;
                od;
                
                #(2) Normalize Pi_mu (all cols have the same norm)
                #Pick the first column of f_mu (since this always exists) and check its norm squared
                f_mu_1 := f_mu{[1..intertwiner_dimension]}{[1]};
                alpha_mu := ComplexNormSquared(f_mu_1); #this is the norm squared of f_nu_1;
                #Divide Pi_mu by the norm squared
                Pi_mu := Pi_mu / alpha_mu ; 
                P := P - Pi_mu;
                #----------------------------Pi_mu Computed and Stored----------------
                

                #--------------------Populate the F_i matrix---------------
                for row_number in [1..intertwiner_dimension] do
                    for col_number in [1..d_c_i] do
                        F_i[row_number][col_number + (mu-1)*d_c_i] := f_mu[row_number][col_number];
                    od;
                od;
                #--------------------------------------------------------
                
                col_index := col_index + d_c_i;
            od;


        fi;
        
        
        #------------Add F_i into the orthogonalized_intertwiner matrix---------
        for row_number in [1..intertwiner_dimension] do
            for col_number in [1..n_i*d_c_i] do
                orthogonalized_intertwiner[row_number][orthogonalized_intertwiner_index + col_number] := F_i[row_number][col_number];
            od;
        od;
        #---------------------------------------------------------------------

        #increment the orthogonalized_intertwiner_index by d_c_i*n_i
        orthogonalized_intertwiner_index := orthogonalized_intertwiner_index + d_c_i*n_i;

    od;

    
    return orthogonalized_intertwiner;

end;

