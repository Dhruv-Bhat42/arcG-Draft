#----------------------------------------------------------------------------------
#Title: A GAP Script to compute the unitary F matrices for the representation category of finite groups 
#Developed by Dhruv Bhat, BSc. Mathematics from New York University Abu Dhabi
#Created under the mentorship of Dr. Sachin Valera, Postdoctoral Associate at the Center for Quantum and Topological Systems, NYU Abu Dhabi
#----------------------------------------------------------------------------------
#Notes to Reader
#This program has been written in GAP, which can be installed from https://www.gap-system.org/install/ 
#In GAP, array indexing starts at 1, not 0, unlike common programming languages like Python. 


#---------LOAD PACKAGES FOR CHAR TABLES, IRREPS AND UNITARISATION COMMANDS-------- 
if not IsPackageLoaded("ctbllib") then
    LoadPackage("ctbllib"); 
fi;

if not IsPackageLoaded("RepnDecomp") then 
    LoadPackage("RepnDecomp");      
fi;

#Loading Dependecies
Read("label_methods.g"); #contains functions to inspect fusion rules and separate multiplicities from irrep labels
Read("randomness.g"); #contains functions to generate random integer matrices
Read("ModifiedGS.g"); #modified Gram Schmidt algorithm 

#WELCOME MESSAGE
Print("Welcome to arcG: Associators of the Representation Categories of finite Groups \n");
Print("This program will output the associators (F symbols) of Rep(G) \nTo begin, please define a finite group \n");

#PROMPTING USER FOR GROUP & RANGE OF RANDOM MATRICES
# G:= 0;
# while not IsGroup(G) do
#     G:= InputFromUser("Your input must be as per GAP syntax: https://docs.gap-system.org/doc/ref/chap50.html \n Sample Input from User: \n SymmetricGroup(5) # this refers to S_5 \n DihedralGroup(10) # this refers to D_5 \n SmallGroup(n,j) #For integers n,j this refers to the jth finite group of order n as per the GAP library \n");
#     if not IsGroup(G) then  
#         Print("Sorry, your input isn't a valid finite group. \n");
#     fi;
# od;

# Print("Thank you! \n Our algorithm makes use of random integer matrices whose entries are sampled uniformly from [-d,d]\\{0} \n A lower value of d would increase compute time, but would allow for a more visually palatable guage for the F symbols. \n");
random_matrix_range := 1;

# while IsInt(random_matrix_range) and (20 < random_matrix_range or random_matrix_range <= 0) do
#     random_matrix_range:= InputFromUser("Please enter your preferred value of d:\n");
#     if not IsInt(random_matrix_range) then  
#         Print("Sorry, your input isn't a valid integer \n");
#     fi;
# od;




#---------DEFINE GROUP----------
#Group Name 
#G := SymmetricGroup(3); #(DEFINE GROUP NAME HERE, OR IN GAP TERMINAL)
Print("\nRunning arcG for Rep(", StructureDescription(G),  "). The GAP I.D. for this group is ", IdGroup(G), "\n \n");

#====================PREPROCESSING BEGINS (Computing and Storing Fusion Rules & Unitary Irreps)==============================
preprocessing_start := Runtime(); #Record GAP runtime at the start of preprocessing
Print("======================Preprocessing =========================\n");

#-----------COMPUTE FUSION RULES---------------------------------------

Read("G_fusion_finder.g");  # This imports fusionRules and irrepNames, and indicates whether multiplicities are present

#------------COMPUTE UNITARISED IRREDUCIBLE REPRESENTATIONS------------------
Read("G_irrep_finder.g"); #This computes unitary irreducible representations for the group, and imports irreps_unitary: a list of all unitary irreps of G sorted in the same order as sortedChars

#----------CREATE A RECORD OF UNITARISED IRREDUCIBLE REPRESENTATIONS AND THEIR RESPECTIVE LABELS FOR EASY ACCESS---------

# Define a "GAP record" (~ dictionary) to store unitary irreps by their labels
irrep_dictionary := rec();
for i in [1..Length(sortedChars)] do
    label := irrepNames[i];  # Get a representation label from the list
    representation := irreps_unitary[i];  # Extract the corresponding representation
    irrep_dictionary.(label) := representation;  # Assign this representation to the record field
od;

#Record GAP runtime at the end of preprocessing
preprocessing_end := Runtime(); 

#Report Preprocessing Runtime
preprocessing_time_taken := preprocessing_end - preprocessing_start;
if preprocessing_end = preprocessing_start then
    Print("[Preprocessing duration: < 1 ms] \n");
else  
    Print("[Preprocessing duration: ", preprocessing_time_taken ," ms] \n");
fi;
#=====================================Preprocessing Complete===========================================

#==================================MAIN ALGORITHM BEGINS==============================================
Print("===========================MAIN ALGORITHM================================ \n");

#-----------(PART 1)- COMPUTING FULL RANK INTERTWINERS [i.e. T(a,b)] FOR EACH NON-TRIVIAL FUSION RULE a*b = ...-----------------------------
intertwiner_computation_start := Runtime();  #Record GAP runtime when the computation of full tank intertwiners begins
Print("Computing full rank intertwiners for each non-trivial fusion rule a*b = ...", "\n");

max_random_matrix_runs := [1, []]; #will store a performance statistic (max num of times random matrices had to be regenerated, and the list of pairs of irreps for which this happened)
T_intertwinerLIST := rec(); #a GAP "record" which works like a python dictionary. 
#The keys will be ordered pairs of irreps (a,b) and the values are the full rank intertwiners T(a,b) that block diagonalize them into a direct sum of irreps


#We iterate over all pairs of irreps and construct unitary intertwiners between their tensor product rep and its block diagonal form
fusionIndex := 1; #to iterate over all fusion rules
for i in [2..repCount] do
    for j in [i..repCount] do
        #We work with irreps_unitary[i] and irreps_unitary[j] in this loop
        identical_irreps := false;
        if i = j then #check if we are fusing identical irreps
            identical_irreps := true;
        fi;

        #-------------CREATE THE BLOCK DECOMPOSITION OF THIS PAIR OF IRREPS---------------
        currentDecomp := fusionRules[fusionIndex]; #a list of fusion outcomes for current pair of irreps
        
        #Create the block diagonal decomposition of the current tensor product representation
        
        #We start by adding the first irrep as many times as its multiplicity
        #We inspect each element of currentDecomp to check the multiplicity and label of the irrep it contains, and add the corresponding irrep (direct sum) to the block diagonal decomp those many times
        
        current_multiplicity := label_info(currentDecomp[1])[1]; #Integer value; multiplicity of the irrep in currentDecomp[k]
        fusion_outcome_label := label_info(currentDecomp[1])[2]; #String- label of the irrep in currentDecomp[k]
        
        block_diagonal_decomp := irrep_dictionary.(fusion_outcome_label); #The block diagonal decomp of the tensor product of the pair of irreps is initialized as whatever the first irrep in the decomposition is. We will add in (direct sum) the others iteratively
        for irrep_copy in [1..current_multiplicity]  do
            if irrep_copy > 1 then #1 copy is already added above, so start adding the second one
                block_diagonal_decomp := DirectSumOfRepresentations([block_diagonal_decomp, irrep_dictionary.(fusion_outcome_label)]);
            fi;
        od;  

        #Add all other irreps
        for k in [1..Length(currentDecomp)] do
            if k > 1 then
                #Split the multiplicity and the label
                current_multiplicity := label_info(currentDecomp[k])[1]; #Integer value; multiplicity of the irrep in currentDecomp[k]
                fusion_outcome_label := label_info(currentDecomp[k])[2]; #String- label of the irrep in currentDecomp[k]
            
                if current_multiplicity = 1 then 
                    #there are no multiplicities. so we add in the current irrep only once
                    block_diagonal_decomp := DirectSumOfRepresentations([block_diagonal_decomp, irrep_dictionary.(currentDecomp[k])]);

                else 
                    #there are multiplicities. so we add this irrep *current_multiplicity* times
                    for irrep_copy in [1..current_multiplicity]  do
                        block_diagonal_decomp := DirectSumOfRepresentations([block_diagonal_decomp, irrep_dictionary.(fusion_outcome_label)]);
                    od;               
                fi;
            fi;
        od;

        #-----APPLY ALGORITHM IN PAPER TO FIND THE UNITARY INTERTWINER T(a,b) THAT RELATES THE TENSOR PRODUCT REP (V_a * V_b) WITH THE BLOCK DIAGONAL REP CONSTRUCTED ABOVE ((+)_c V_c)---
        tensor_rep_dimension := DegreeOfRepresentation(irreps_unitary[i])*DegreeOfRepresentation(irreps_unitary[j]); #the dimension of the tensor product representation in this iteration of the loop
        
        #(1)- Find unitary intertwiner that block diagonalizes irreps_unitary[i]*irreps_unitary[j]
        sakata1 := NullMat(tensor_rep_dimension, tensor_rep_dimension);

        
        #Since we use a random matrix in this algorithm, we may need to change the choice of random matrix if sakata1 happens to be singular
        iteration_count := 0; #a safety measure to prevent an infinite loop 
        
        while IsZero(DeterminantMat(sakata1)) do
            randMatrix1 := GenerateRandomMatrix(tensor_rep_dimension, -1*random_matrix_range, random_matrix_range);  # Random square matrix of dimension matching the tensor product representation populated by integers between the specified range
            for g in Elements(G) do
                sakata1:= sakata1 + KroneckerProduct(irreps_unitary[i](Inverse(g)), irreps_unitary[j](Inverse(g)))*randMatrix1*block_diagonal_decomp(g);
            od;
            iteration_count := iteration_count + 1;
            
            if iteration_count > 20 and IsZero(DeterminantMat(sakata1)) then
                Print("WARNING: Intertwiner for ", irrepNames[i], " tensor ", irrepNames[j],  " was not full rank for 20 runs! \n(Script is potentially stuck in a loop, since this step relies on random matrices) \n");    
            fi;

            #---Update the statistics on the maximum number of runs needed for this loop over all pairs of intertwiners
            if iteration_count > max_random_matrix_runs[1] then #if higher number of max iterations found
                max_random_matrix_runs[1] := iteration_count; #update the max_random_matrix_runs variable
                max_random_matrix_runs[2] := []; #start a new empty list of irrep pairs that require these many re-runs                
            fi;
            if iteration_count = max_random_matrix_runs[1] then #if this iteration count matches max_random_matrix_runs thus far
                Add(max_random_matrix_runs[2], Concatenation(irrepNames[i],"*" ,irrepNames[j])); #add this irrep pair to max_random_matrix_runs[2] (i.e. the list of irrep pairs which need the most amount of re-runs)           
            fi;
            #---Updating Complete

        od;
        
        #Label and store this matrix
        intertwiner_label1 := Concatenation("T_",irrepNames[i], "*", irrepNames[j]); #the label for the intertwiner that block diagonalizes irrepNames[i] \tensor irrepNames[j]
        T_intertwinerLIST.(intertwiner_label1) := sakata1; #add this intertwiner to the record (dictionary) of intertwiners with the appropriate label
        

        #(2)- Find unitary intertwiner that block diagonalizes irreps_unitary[j]*irreps_unitary[i]
        if not identical_irreps then #only if the two irreps are distinct
            sakata2 := NullMat(tensor_rep_dimension, tensor_rep_dimension);
            iteration_count := 0; #a safety measure to prevent an infinite loop
            
            while IsZero(DeterminantMat(sakata2)) do
                randMatrix2 := GenerateRandomMatrix(tensor_rep_dimension, -1*random_matrix_range, random_matrix_range);  # Random invertible square matrix of dimension matching the tensor product representation populated by integers
                for g in Elements(G) do
                    sakata2:= sakata2 + KroneckerProduct(irreps_unitary[j](Inverse(g)), irreps_unitary[i](Inverse(g)))*randMatrix2*block_diagonal_decomp(g);
                od;
                iteration_count := iteration_count + 1;
                
                if iteration_count >= 20 and IsZero(DeterminantMat(sakata2)) then
                    Print("WARNING: Intertwiner for ", irrepNames[j], " tensor ", irrepNames[i],  " was not full rank for 20 runs! \n(Script is potentially stuck in a loop, since this step relies on random matrices) \n");
                fi;

                #---Update the statistics on the maximum number of runs needed for this loop over all pairs of intertwiners
                if iteration_count > max_random_matrix_runs[1] then #if higher number of max iterations found
                    max_random_matrix_runs[1] := iteration_count; #update the max_random_matrix_runs variable
                    max_random_matrix_runs[2] := []; #start a new empty list of irrep pairs that require these many re-runs                
                fi;
                if iteration_count = max_random_matrix_runs[1] then #if this iteration count matches max_random_matrix_runs thus far
                    Add(max_random_matrix_runs[2], Concatenation(irrepNames[j],"*" ,irrepNames[i])); #add this irrep pair to max_random_matrix_runs[2] (i.e. the list of irrep pairs which need the most amount of re-runs)           
                fi;
                #---Updating Complete   
            od;
            
            #Label and Store this matrix
            intertwiner_label2 := Concatenation("T_", irrepNames[j], "*", irrepNames[i]); #the dictionary key for the intertwiner that block diagonalizes irrepNames[j] \tensor irrepNames[i]
            T_intertwinerLIST.(intertwiner_label2) := sakata2; #add this intertwiner to the record (dictionary) of intertwiners with the appropriate label


        fi;

        #increment the index to work with the next fusion rule
        fusionIndex := fusionIndex + 1; 
    od;
od;

#Record GAP runtime when the computation of full tank intertwiners ends
intertwiner_computation_end := Runtime(); 

#Print success message & total runtime for the computation of full rank intertwiners
intertwiner_computation_time := intertwiner_computation_end - intertwiner_computation_start;
if intertwiner_computation_time = 0 then 
    Print("All full rank intertwiners computed! [< 1 ms] \n");
else     
    Print("All full rank intertwiners computed! [", intertwiner_computation_time ," ms] \n");
fi;


#---------(Part 2) - ORTHOGONALISING COLUMNS OF ALL FULL RANK INTERTWINERS-------------------
Print("Orthogonalising columns of full rank interwiners...\n");
orthogonalization_start_time := Runtime(); #Record GAP runtime when the column orthogonalization of full tank intertwiners begins
u_intertwinerLIST := rec();
#The keys will be ordered pairs of irreps and the values are the unitary intertwiners that block diagonalize them into a direct sum of irreps

fusionIndex := 1;
for i in [2..repCount] do
    for j in [i..repCount] do
        
        #We work with irreps_unitary[i] and irreps_unitary[j] in this loop
        identical_irreps := false;
        if i = j then #check if we are fusing identical irreps
            identical_irreps := true;
        fi;
        
        #Get a list of fusion outcomes for current pair of irreps
        currentDecomp := fusionRules[fusionIndex]; 


        intertwiner1 := Concatenation("T_", irrepNames[i], "*", irrepNames[j]); #the dictionary key for the intertwiner that block diagonalizes irrepNames[i] \tensor irrepNames[j]
        intertwiner1_matrix := T_intertwinerLIST.(intertwiner1) ; 

        #Apply Modified Gram Schmidt to intertwiner1
        intertwiner1_matrix := Modified_Gram_Schmidt(intertwiner1_matrix, currentDecomp);
        #VALIDATION CHECK FOR ORTHOGONALITY
        #Print("Are Cols Orthogonal?: ", AreColumnsOrthogonal(sakata1), "\n");
        
        unitarized_intertwiner1_label := Concatenation("U_", irrepNames[i], "*", irrepNames[j]);
        #Add the entry in the dictionary of unitary intertwiners
        u_intertwinerLIST.(unitarized_intertwiner1_label) :=  intertwiner1_matrix; 
        
        #Print(unitarized_intertwiner1_label, ": ", intertwiner1_matrix, "\n");

        if not identical_irreps then
            intertwiner2 := Concatenation("T_", irrepNames[j], "*", irrepNames[i]); #the dictionary key for the intertwiner that block diagonalizes irrepNames[j] \tensor irrepNames[i]
            intertwiner2_matrix := T_intertwinerLIST.(intertwiner2);

            #Apply Modified Gram Schmidt to intertwiner2
            intertwiner2_matrix := Modified_Gram_Schmidt(intertwiner2_matrix, currentDecomp);
            #VALIDATION CHECK FOR ORTHOGONALITY
            #Print("Are Cols Orthogonal?: ", AreColumnsOrthogonal(sakata1), "\n");

            unitarized_intertwiner2_label :=  Concatenation("U_", irrepNames[j], "*", irrepNames[i]);
            #Add the entry in the dictionary of unitarized intertwiners
            u_intertwinerLIST.(unitarized_intertwiner2_label) :=  intertwiner2_matrix; 
            #Print(unitarized_intertwiner2_label, ": ", intertwiner2_matrix, "\n");
        fi;

        fusionIndex := fusionIndex + 1;
    od;
od;

orthogonalization_end_time := Runtime(); #Record GAP runtime when the column orthogonalization of full tank intertwiners begins

#Report success message & runtime for the column orthogonalization of full tank intertwiners
orthogonalization_time_taken := orthogonalization_end_time - orthogonalization_start_time;
if orthogonalization_time_taken = 0 then
    Print("Orthogonalization Complete! [< 1 ms] \n");
else  
    Print("Orthogonalization Complete! [", orthogonalization_time_taken," ms] \n");
fi;

#-------------------------------EXTRACT THE TRIVALENT VERTICES IN REP(G)-----------------------------
#This part of the code slices portions of each of the "Unitary Intertwiners" and assigns them to trivalent vertices
trivalent_vertex_extraction_start := Runtime(); #Record GAP runtime when extraction of trivalent vertices begins
Print("Commencing Trivalent Vertex Extraction! \n");
#The trivalent vertices dictionary (splitting vertices)
trivalent_vertices := rec();
#Key-Value Pairs are of 2 types
#Type 1 Keys are of form a*b_c__mu and the values are the matrix form of the splitting vertex associated to the key
#for example if R_2 * R_3_1 = 2R_3_2 then
#then R_2*R_3_1_R_3_2__1 would be the key for the matrix associated to the linear map represented by the splitting vertex:
#      | R_3_2
#      | 
#      | 1
# R_2 / \ R_3_1
#Type 2 Keys are of form NORMSQ_a*b_c__mu and the corresponding value is a scalar value, which is the common norm-squared of all the columns in the matrix associated to key a*b_c__mu

fusionIndex := 1;

for i in [2..repCount] do
    for j in [i..repCount] do
        #We will define all (splitting) trivalent vertices pertaining to irreps_unitary[i] and irreps_unitary[j]
        #-----
        currentDecomp := fusionRules[fusionIndex]; #The list of the fusion outcomes of irreps_unitary[i] tensor irreps_unitary[j]
        #-----
        #Check if we are fusing identical irreps
        identical_irreps := false;
        if i = j then 
            identical_irreps := true;
        fi;
        #----

        #(1)- GET ALL (splitting) TRIVALENT VERTICES for the ORDERED PAIR (irreps_unitary[i],irreps_unitary[j])
        current_unitary_intertwiner := u_intertwinerLIST.(Concatenation("U_", irrepNames[i],"*",irrepNames[j])); #access the unitary intertwiner that block diagonalizes irreps_unitary[i] \tensor irreps_unitary[j] into a direct sum of irreps
        
        u_size := Length(current_unitary_intertwiner); #this is the dimension of the current_unitary_intertwiner matrix
    
        #Create a variable to record the number of columns of current_unitary_intertwiner which have been extracted
        u_col_index := 1; #each time a trivalent vertex is extracted, this variable will be incremented by the dimension of the fusion outcome (the dimension of the irrep)

        #Iterate over each fusion outcome & multiplicity to create each trivalent vertex
        for k in [1..Length(currentDecomp)] do
            
            #Record the label, multiplicity and degree of the fusion outcome in this iteration
            current_multiplicity := label_info(currentDecomp[k])[1]; #Integer value; multiplicity of the irrep in currentDecomp[k]
            fusion_outcome_label := label_info(currentDecomp[k])[2]; #String- label of the irrep in currentDecomp[k]
            irrep_dimension := DegreeOfRepresentation(irrep_dictionary.(fusion_outcome_label)); #degree of irrep contained in currentDecomp[k]

            # if current_multiplicity = 1 then
            #     Print("Extracting the trivalent vertex for the fusion ", irrepNames[i], "*", irrepNames[j], " into ", fusion_outcome_label, "...", "\n");
            # else
            #     Print("Extracting the indexed trivalent vertices for the fusion ", irrepNames[i], "*", irrepNames[j], " into ", fusion_outcome_label, "...", "\n");
            # fi;

            #We have to extract *current_multiplicity* trivalent vertices for the splitting of fusion_outcome_label into irreps_unitary[i] & irreps_unitary[j]
            for multiplicity_index in [1..current_multiplicity] do #Loop over the multiplicity
                
                #Create a label for the trivalent vertex associated to this indexed triplet of irreps in the form "irrepNames[i]*irrepNames[j]_(fusion_outcome_label)__multiplicity"
                trivalent_vertex_label := Concatenation(irrepNames[i], "*", irrepNames[j], "_", fusion_outcome_label, "__", String(multiplicity_index));

                #Pick the appropriate columns of the unitary intertwiner that block diagonalizes irrepNames[i] and irrepNames[j] 
                if irrep_dimension = 1 then
                    vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index]}; #this extracts only the column at index u_col_index (Note: GAP indexing starts at 1)
                else
                    vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index..(u_col_index+irrep_dimension-1)]}; #this extracts all columns from col u_col_index and (u_col_index+irrep_dimension);
                fi;
                
                #Associate this vertex label with the vertex data in the trivalent vertex record (which works like a dictionary with key/value pairs)
                trivalent_vertices.(trivalent_vertex_label) := vertex_data;
                
                #Add the common norm squared of the cols of vertex_data into the dictionary as well
                vertex_norm_squared := ComplexNormSquared(vertex_data{[1..u_size]}{[1]}); #norm of the first col of vertex_data
                vertex_norm_squared_label := Concatenation("NORMSQ_", trivalent_vertex_label);
                trivalent_vertices.(vertex_norm_squared_label) := vertex_norm_squared;

                #DEBUGGING
                # Print(trivalent_vertex_label, ": \n");
                # Print(vertex_data, "\n");
                # Print(vertex_norm_squared_label, ": ", vertex_norm_squared, "\n");
                u_col_index := u_col_index + irrep_dimension ; #check if col_index was correct
                
            od;
        od;
        

        #(2)- GET ALL (splitting) TRIVALENT VERTICES for the ORDERED PAIR (irreps_unitary[j],irreps_unitary[i])
        if not identical_irreps then #as long as we have distinct irreps
      
            current_unitary_intertwiner := u_intertwinerLIST.(Concatenation("U_", irrepNames[j],"*",irrepNames[i])); ##access the unitary intertwiner that block diagonalizes irreps_unitary[j] \tensor irreps_unitary[i] into a direct sum of irreps
            
            #Create a variable to record the number of columns of current_unitary_intertwiner which have been extracted
            u_col_index := 1; #each time a trivalent vertex is extracted, this variable will be incremented by the dimension of the fusion outcome (the dimension of the irrep)

            for k in [1..Length(currentDecomp)] do
                #Record the label, multiplicity and degree of the fusion outcome in this iteration
                current_multiplicity := label_info(currentDecomp[k])[1]; #Integer value; multiplicity of the irrep in currentDecomp[k]
                fusion_outcome_label := label_info(currentDecomp[k])[2]; #String- label of the irrep in currentDecomp[k]
                irrep_dimension := DegreeOfRepresentation(irrep_dictionary.(fusion_outcome_label)); #degree of irrep contained in currentDecomp[k]

                #DEBUGGING
                # if current_multiplicity = 1 then
                #     Print("Extracting the trivalent vertex for the fusion ", irrepNames[j], "*", irrepNames[i], " into ", fusion_outcome_label, "...", "\n");
                # else
                #     Print("Extracting the indexed trivalent vertices for the fusion ", irrepNames[j], "*", irrepNames[i], " into ", fusion_outcome_label, "...", "\n");
                # fi;

                #We have to extract *current_multiplicity* trivalent vertices for the splitting of fusion_outcome_label into irreps_unitary[i] & irreps_unitary[j]
                for multiplicity_index in [1..current_multiplicity] do #Loop over the multiplicity
                    
                    #Create a label for the trivalent vertex associated to this indexed triplet of irreps in the form "irrepNames[i]*irrepNames[j]_(fusion_outcome_label)__multiplicity"
                    trivalent_vertex_label := Concatenation(irrepNames[j], "*", irrepNames[i], "_", fusion_outcome_label, "__", String(multiplicity_index));

                    #Pick the appropriate columns of the unitary intertwiner that block diagonalizes irrepNames[j] tensor irrepNames[i] 
                    if irrep_dimension = 1 then
                        vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index]}; #this extracts only the column at index u_col_index (Note: GAP indexing starts at 1)
                    else
                        vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index..(u_col_index+irrep_dimension-1)]}; #this extracts all columns from col u_col_index and (u_col_index+irrep_dimension);
                    fi;
                    
                    #Associate this vertex label with the vertex data in the trivalent vertex record (which works like a dictionary with key/value pairs)
                    trivalent_vertices.(trivalent_vertex_label) := vertex_data;
                    
                    #Add the common norm squared of the cols of vertex_data into the dictionary as well
                    vertex_norm_squared := ComplexNormSquared(vertex_data{[1..u_size]}{[1]}); #norm of the first col of vertex_data
                    vertex_norm_squared_label := Concatenation("NORMSQ_", trivalent_vertex_label);
                    trivalent_vertices.(vertex_norm_squared_label) := vertex_norm_squared;

                    u_col_index := u_col_index + irrep_dimension ;
                    
                od;
            od;
            
        fi;
        
        #Increment fusionIndex to work with the next fusion rule
        fusionIndex := fusionIndex + 1;

    od;
od;

#Also add all the trivial trivalent vertices (fusion of trivial irrep with irrep x into irrep x itself)
for i in [1..repCount] do
        degree_i := DegreeOfRepresentation(irreps_unitary[i]);
        trivial_label1 := Concatenation("R_I", "*", irrepNames[i], "_", irrepNames[i], "__1"); #Fusion of R_I with irrep i into irrep i
        trivial_label2 := Concatenation(irrepNames[i], "*", "R_I", "_", irrepNames[i], "__1"); #Fusion of irrep i with R_I into irrep i
        #Both the trivial trivalent vertices are simply the identity map on the representation space of irrep i
        trivalent_vertices.(trivial_label1) := IdentityMat(degree_i); 
        trivalent_vertices.(trivial_label2) := IdentityMat(degree_i);

        #Add the common norm squared of the cols of vertex_data into the dictionary as well
        trivalent_vertices.(Concatenation("NORMSQ_", trivial_label1)) := 1; 
        trivalent_vertices.(Concatenation("NORMSQ_", trivial_label2)) := 1;


od;

#Record GAP runtime when extraction of trivalent vertices begins
trivalent_vertex_extraction_end := Runtime();

#Print success message & time taken for extraction of trivalent vertices
trivalent_vertex_extraction_time := trivalent_vertex_extraction_end - trivalent_vertex_extraction_start;
if trivalent_vertex_extraction_time = 0 then   
    Print("Trivalent Vertex Extraction Complete! [< 1 ms] \n");
else
    Print("Trivalent Vertex Extraction Complete! [", trivalent_vertex_extraction_time ," ms] \n");
fi;

#------------------------------------COMPUTE F SYMBOLS VIA THE GRAPHICAL CALCULUS OF A FUSION CATEGORY---------------------------------
#-----Create a record (dictionary) of fusion rules for on-command access to fusion rules (instead of iterating through every rule) ----
fusionRules_record := rec(); #dictionary of fusion rules in Rep(G). 
#The Key-Value pairs in fusionRules_record are of 2 Types (for cleaner data extraction later)
#A Type 1 Key is of form "irrep a*irrep b__irrep c", where c is a valid fusion outcome. The value is the multiplicity of c (i.e. N^{ab}_c)
#A Type 2 Key is of form "irrep a*irrep b". The value is a list of the labels of fusion outcomes of a and b (with no multiplicity information)

#We will use the fusionRules list to construct this dictionary.
fusionIndex := 1; #to access each element of fusionRules sequentially and populate the dictionary 

for i in [2..repCount] do
    for j in [i..repCount] do
        fusion_outcome_list := []; #A list of fusion outcomes of irreps i and j (without multiplicities)
        currentDecomp := fusionRules[fusionIndex]; #currentDecomp is a list of fusion outcomes of i and j containing multiplicities as well
        
        #Extract distinct fusion outcomes and their respective multiplicities
        for outcome in currentDecomp do 
            fusion_outcome_multiplicity := label_info(outcome)[1]; #this is the multiplicity of the fusion outcome
            fusion_outcome_label := label_info(outcome)[2]; #this is the label of the fusion outcome
            
            #Add the label of this fusion outcome to fusion_outcome_list
            Add(fusion_outcome_list, fusion_outcome_label); 

            #---------Add Key-Value Pairs of Type 1----------
            #Keys of Type 1 are pairs of Irreps. Their values will be distinct fusion outcomes
            multiplicity_key_1 := Concatenation(irrepNames[i], "*", irrepNames[j], "__", fusion_outcome_label); #key for i*j__(fusion_outcome_label)
            multiplicity_key_2 := Concatenation(irrepNames[j], "*", irrepNames[i], "__", fusion_outcome_label); #key for i*j__(fusion_outcome_label)
            #Assign their values (Multiplicity of fusion_outcome_label in i*j and j*i)
            fusionRules_record.(multiplicity_key_1) := fusion_outcome_multiplicity;
            fusionRules_record.(multiplicity_key_2) := fusion_outcome_multiplicity;
        od;

    	#---------Add Key-Value Pairs of Type 2----------
        #Keys of Type 2 are pairs of Irreps. Their values will be distinct fusion outcomes
        fusion_key_1 := Concatenation(irrepNames[i], "*", irrepNames[j]); #key for i*j
        fusion_key_2 := Concatenation(irrepNames[j], "*", irrepNames[i]); #key for j*i
        #Assign their values (fusion rules for i*j and j*i are always the same- this is a commutative fusion ring)
        fusionRules_record.(fusion_key_1) := fusion_outcome_list;
        fusionRules_record.(fusion_key_2) := fusion_outcome_list;

        #Increment the index to work with the next pair of irreps (in the next iteration of the loop)
        fusionIndex := fusionIndex + 1;
    od;
od;

#Also add trivial fusion rules into fusionRules_record for convenience
for i in [2..repCount] do
    #Key-Value Pairs of Type 1
    trivial_multiplicity_label1 := Concatenation("R_I", "*", irrepNames[i], "__", irrepNames[i]);
    trivial_multiplicity_label2 := Concatenation(irrepNames[i] ,"*", "R_I", "__", irrepNames[i]);
    fusionRules_record.(trivial_multiplicity_label1):= 1; #Since N^{Ix}_x = 1 for all x
    fusionRules_record.(trivial_multiplicity_label2):= 1; #Since N^{xI}_x = 1 for all x
    
    #Key-Value Pairs of Type 2
    trivial_rule_label1 := Concatenation("R_I", "*", irrepNames[i]);
    trivial_rule_label2 := Concatenation(irrepNames[i] ,"*", "R_I");
    fusionRules_record.(trivial_rule_label1):= [irrepNames[i]]; #Since (R_I * irrep i = irrep i) is the entirety of the fusion rule
    fusionRules_record.(trivial_rule_label2):= [irrepNames[i]]; #Since (irrep i * R_I = irrep i) is the entirety of the fusion rule
od;
#-------------------------------------------------------

#----Compute the Fusion Outcomes of Every Distinct Triplet of Irreps of the Group G----
triplet_fusion_outcomes := rec(); #a dictionary of fusion outcomes for every ordered triplet of irreps in Rep(G). 
#Each key is an ordered tuple of non-trivial irreps and the corresponding value is a list of distinct fusion outcomes of the triplet

fusionIndex := 1;
for a in [2..repCount] do #first irrep label a
    for b in [a..repCount] do #second irrep label b
        #Obtain a list of distinct fusion outcomes of a*b and store it in currentDecomp
        currentDecomp := [];
        #Use the fusionRules_record dictionary to populate currentDecomp with the distinct fusion outcomes of a*b
        fusion_key := Concatenation(irrepNames[a], "*", irrepNames[b]); #This key can access the fusion rules for a*b via the fusionRules_record dictionary
        currentDecomp := fusionRules_record.(fusion_key);
        #------
        for c in [2..repCount] do #third irrep label c (iterate over all other non-trivial irreps)
            #---We will now find the list of all distinct fusion outcomes of a*b*c----
            abc_distinct_outcomes := []; #initialize a list to store all distinct outcomes
            #------
            for k in [1..Length(currentDecomp)] do #iterate over all fusion outcomes of a and b
                ab_fusion_outcome := currentDecomp[k]; #store the specific fusion outcome of a and b listed at currentDecomp[k]
                
                #Find all the fusion outcomes of currentDecomp[k] with irrep c
                ab_k_c_outcomeList := []; #a list of all fusion outcomes of currentDecomp[k] and c. Note on variable name: "ab_k" reminds us that currentDecomp[k] was a fusion outcome of a & b
                #Use the fusionRules_record dictionary to populate ab_k_c_outcomeList
                k_c_outcomeList_label := Concatenation(ab_fusion_outcome, "*", irrepNames[c]); #This key can access fusion rules of currentDecomp[k] and c via the fusionRules_record dictionary               
                ab_k_c_outcomeList := fusionRules_record.(k_c_outcomeList_label); 
                
                #Each fusion outcome of currentDecomp[k] and irrep c is a fusion outcome of a*b*c; Thus, we populate abc_distinct_outcomes via ab_k_c_outcomeList
                for outcome in ab_k_c_outcomeList do
                    if not outcome in abc_distinct_outcomes then #check if this is a new a*b*c fusion outcome that has not been recorded before
                        Add(abc_distinct_outcomes, outcome);
                    fi;
                od;
            od;
            
            #STORE THE DATA IN THE triplet_fusion_outcomes record (all orders of the triplet {a,b,c} are covered through this implementation)

            triplet_label1 := Concatenation(irrepNames[a], "*", irrepNames[b], "*", irrepNames[c]);
            triplet_fusion_outcomes.(triplet_label1) := abc_distinct_outcomes;
            #Print("The outcomes of ", triplet_label1, " are: ", abc_distinct_outcomes, "\n");
            
            #If irrep a and b are distinct, then store the fusion rule for the swapped pair (b,a) to cover all permutations through which this fusion rule can be accessed
            if not a = b then 
                triplet_label2 := Concatenation(irrepNames[b], "*", irrepNames[a], "*", irrepNames[c]);
                triplet_fusion_outcomes.(triplet_label2) := abc_distinct_outcomes;
                #Print("The outcomes of ", triplet_label2, " are: ", abc_distinct_outcomes, "\n");
                
            fi;
        od; 

        fusionIndex := fusionIndex + 1; #Increment to work with the next fusion rule
    od;
od;
#--------------------------------------------------

#--------Label and Compute F symbols--------------
Fsymbol_computation_start := Runtime(); #Record GAP runtime when computation of F symbols begins
Print("Computing the F symbols of Rep(G) through the graphical calculus!...", "\n");
F_symbol_record := rec(); #this is the dictionary of F symbols of Rep(G)
#A key in this record is the label of an F symbol. This label will be a list of 6 labels {a,b,c,e,f,d} which refers to entry (f,d) of the F-matrix [F^{abc}_e] (fusion of irreps a, b and c into e)
#Naturally, the value of each key will be the numerical value of the F symbol. 

# We will now methodically LABEL and COMPUTE every non-trivial F symbol of Rep(G)
#No need to consider F symbols involving the trivial irrep- this is set to 1 by definition in a fusion category! Thus we start indexing at 2. 

#Access Appropriate Filename & File Stream for Output Data
filename := Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Output_Data.txt");
stream := OutputTextFile(filename, true); #true - append mode on. false would mean create new file/empty existing one
#We must use append mode here, since G_fusion_finder.g has already printed in the fusion rules into this file
SetPrintFormattingStatus( stream, false ); #makes sure there are no arbitrary line breaks
PrintTo(stream, "\nHere are the F Symbols for ", StructureDescription(G), "\n");

total_non_trivial_F_matrices := 0; #will store the number of non-trivial F matrices computed (i.e. F^{abc}_d wherein a,b and c are not the identity)
total_non_trivial_F_symbols := 0; #will store the total number of non-trivial F symbols computed

for a in [2..repCount] do #pick first irrep label irrepNames[a]
    for b in [2..repCount] do #pick second irrep label irrepNames[b]
        for c in [2..repCount] do #pick third irrep label irrepNames[c] (here iterate over all non-trivial irreps)
            #We are working with V^{abc} in this iteration of the loop
            #We want to compute F matrices for the fusion spaces V^{abc}_d for all valid d

            #First we find all valid d (i.e. access all the valid fusion outcomes of a,b and c via the triplet_fusion_record dictionary)
            triplet_label := Concatenation(irrepNames[a], "*", irrepNames[b], "*", irrepNames[c]);
            abc_distinct_outcomes := triplet_fusion_outcomes.(triplet_label); #access all the distinct fusion outcomes of a,b and c 
            
            #Now, we iterate over all distinct fusion outcomes of a,b and c V^{abc}_d for all d
            for abc_outcome in abc_distinct_outcomes do
                #Now we are ready to work with V^{abc}_d for a particular, valid choice of d (note: d is represented by the "abc_outcome" variable in the current iteration of the loop) 
                total_non_trivial_F_matrices := total_non_trivial_F_matrices + 1;
                #We will compute labels and entries of the F matrix F^{abc}_d, but for this we must first find the left and right hand basis of V^{abc}_d
                #Print("Extracting the entries of the F matrix F^{", triplet_label, "}_", abc_outcome, "\n");
            
                #(1)- FIND THE LEFT AND RIGHT HAND FUSION BASIS FOR V^{abc}_d
                leftBasis := []; #will store fusion outcomes of a and b that can produce d upon fusion with c (d is the value of abc_outcome in this loop)
                rightBasis := []; #will store fusion outcomes of b and c that can produce d upon fusion with a (d is the value of abc_outcome in this loop)

                #COMPUTING LEFT BASIS- Fuse a and b first
                ab_outcomes := fusionRules_record.(Concatenation(irrepNames[a], "*", irrepNames[b])); #access the fusion outcomes of a*b
                for ab_outcome in ab_outcomes do #for an outcome e of a*b (i.e. irrep e is the iteration variable ab_outcome in this loop)
                    e_c_label := Concatenation(ab_outcome, "*", irrepNames[c]); #the dictionary key for the list of fusion outcomes of e*c
                    ec_outcomeList := fusionRules_record.(e_c_label); #access the fusion outcomes of irrep d and c
                    if abc_outcome in ec_outcomeList then #if d (abc_outcome in this loop) appears amongst the fusion outcomes of e*c (where e is ab_outcome in this loop)
                        Add(leftBasis, ab_outcome); #add e (ab_outcome in this loop) to leftBasis (left hand basis list of V^{abc}_d)
                    fi;
                od;

                #COMPUTING RIGHT BASIS- Fuse b and c first
                bc_outcomes := fusionRules_record.(Concatenation(irrepNames[b], "*", irrepNames[c])); #access the fusion outcomes of b*c
                for bc_outcome in bc_outcomes do #for an outcome f of b*c (f is the iteration variable bc_outcome in this loop)
                    a_f_label := Concatenation(irrepNames[a], "*", bc_outcome); #the dictionary key for the list of the fusion outcomes of a*f
                    af_outcomeList := fusionRules_record.(a_f_label); #access the list of fusion outcomes of irrep a and f
                    if abc_outcome in af_outcomeList then #if d (abc_outcome in this loop) appears amongst the fusion outcomes of a*f (where f is bc_outcome in this loop)
                        Add(rightBasis, bc_outcome); #add f (bc_outcome in this loop) to rightBasis (right hand basis list of V^{abc}_e)
                    fi;
                od;

                #(2)---------------Compute LEFT and RIGHT HAND TREES in V^{abc}_d ------------------------- 
                #Streamline variable names
                irrep_a :=  irrepNames[a];
                irrep_b := irrepNames[b];
                irrep_c := irrepNames[c];
                irrep_d := abc_outcome;
                
                #----LEFT HAND TREES----
                #We first compute all binary trees of the following form: (DIAGRAM 1)
                #         |
                #  (nu_L) d
                #        / \
                #(mu_L) e   \
                #      / \   \
                #     a   b   c
                #Note: Time flows bottom to top. If morphisms f1,f2 are represented by diagrams d1, d2 respectively then stacking d2 on d1 represents morhpism "f2 of f1"
                leftTrees_record := rec(); #A dictionary of all left hand trees. Keys are 3-tuples (e, mu, nu). Values are the matrix form of linear operator this diagram represents
                #The keys for the dictionary are of form "e,mu_L,nu_L" where mu_L and nu_L in DIAGRAM 1 will take on values {1...N^{ab}_e} and {1...N^{ec}_d} respectively
                leftTrees_ordered_keys := []; #An array containing all keys of leftTrees_record in the appropriate order (Convention 4.1) 
                

                for irrep_e in leftBasis do
                    #Find N^{ab}_e and N^{ec}_d to access the relevant trivalent vertices
                    #We use key-value pairs of Type 1 in the fusionRules_record dictionary
                    ab_e_multiplicity_key := Concatenation(irrep_a, "*", irrep_b, "__", irrep_e); 
                    N_ab_e := fusionRules_record.(ab_e_multiplicity_key); #N^{ab}_e
                    ec_d_multiplicity_key := Concatenation(irrep_e, "*", irrep_c, "__", irrep_d); 
                    N_ec_d := fusionRules_record.(ec_d_multiplicity_key); #N^{ec}_d
                    
                    #Extract the relevant indexed trivalent vertices and create the dictionary record for each left hand tree in V^{abc}_d
                    for mu_L in [1..N_ab_e] do
                        for nu_L in [1..N_ec_d] do
                            #Create the label of this left hand tree
                            leftTree_label := Concatenation(irrep_e, ",", String(mu_L), ",", String(nu_L));
                            leftTree_norm_squared_label := Concatenation("NORMSQ_", irrep_e, ",", String(mu_L), ",", String(nu_L));
                            #Compute the matrix form of the linear operator it represents
                            c_degree := DegreeOfRepresentation(irreps_unitary[c]); #dimension of irrep c
                            c_identity := IdentityMat(c_degree); #identity map on the representation space of irrep c

                            ab_e_mu_L_fusion_vertex_label := Concatenation(irrep_a, "*", irrep_b, "_", irrep_e, "__", String(mu_L));
                            ab_e_mu_L_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(ab_e_mu_L_fusion_vertex_label))); #Hermitian Conjugate of the Splitting Vertex e -> ab through mu_L
                            #Print("ab_e_mu_L fusion vertex is : ", ab_e_mu_L_fusion_vertex, "\n");
                            ec_d_nu_L_fusion_vertex_label := Concatenation(irrep_e, "*", irrep_c, "_", irrep_d, "__", String(nu_L));
                            ec_d_nu_L_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(ec_d_nu_L_fusion_vertex_label))); #Hermitian Conjugate of the Splitting Vertex d -> ec through nu_L
                            #Compute diagram 1 via diagram stacking 
                            diagram1 := ec_d_nu_L_fusion_vertex*KroneckerProduct(ab_e_mu_L_fusion_vertex, c_identity);
                            diagram1_norm_squared := [trivalent_vertices.(Concatenation("NORMSQ_", ab_e_mu_L_fusion_vertex_label)) , trivalent_vertices.(Concatenation("NORMSQ_", ec_d_nu_L_fusion_vertex_label))];
                            #Add this left hand tree into the dictionary
                            leftTrees_record.(leftTree_label) := diagram1;
                            leftTrees_record.(leftTree_norm_squared_label) := diagram1_norm_squared;

                            #Add the key corresponding to this tree in the leftTrees_ordered_keys array
                            Add(leftTrees_ordered_keys, leftTree_label);

                        od;
                    od;
                od;

                #----RIGHT HAND TREES----
                #We first compute all binary trees of the following form: (DIAGRAM 2)
                #     |
                #     d (nu_R)
                #    / \
                #   /   f (mu_R)
                #  /   / \
                # a   b   c   
                #Note: Time flows bottom to top. If morphisms f1,f2 are represented by diagrams d1, d2 respectively then stacking d2 on d1 represents morhpism "f2 of f1"
                rightTrees_record := rec(); #A dictionary of all left hand trees. Keys are 3-tuples (e, mu, nu). Values are the matrix form of linear operator this diagram represents
                #The keys for the dictionary are of form "e,mu_R,nu_R" where mu_R and nu_R in DIAGRAM 1 will take on values {1...N^{bc}_f} and {1...N^{af}_d} respectively
                
                rightTrees_ordered_keys := []; #An array containing all keys of rightTrees_record in the appropriate order (Convention 4.1) 
                
                for irrep_f in rightBasis do
                    #Find N^{bc}_f and N^{af}_d to access the relevant trivalent vertices
                    #We use key-value pairs of Type 1 in the fusionRules_record dictionary
                    bc_f_multiplicity_key := Concatenation(irrep_b, "*", irrep_c, "__", irrep_f); 
                    N_bc_f := fusionRules_record.(bc_f_multiplicity_key); #N^{bc}_f
                    af_d_multiplicity_key := Concatenation(irrep_a, "*", irrep_f, "__", irrep_d); 
                    N_af_d := fusionRules_record.(af_d_multiplicity_key); #N^{af}_d
                    
                    #Extract the relevant indexed trivalent vertices and create the dictionary record for each left hand tree in V^{abc}_d
                    for mu_R in [1..N_bc_f] do
                        for nu_R in [1..N_af_d] do
                            #Create the label of this left hand tree
                            rightTree_label := Concatenation(irrep_f, ",", String(mu_R), ",", String(nu_R));
                            rightTree_norm_squared_label := Concatenation("NORMSQ_", irrep_f, ",", String(mu_R), ",", String(nu_R));
                            #Compute the matrix form of the linear operator it represents
                            a_degree := DegreeOfRepresentation(irreps_unitary[a]); #dimension of irrep a
                            a_identity := IdentityMat(a_degree); #identity map on the representation space of irrep a

                            bc_f_mu_R_fusion_vertex_label:= Concatenation(irrep_b, "*", irrep_c, "_", irrep_f, "__", String(mu_R));
                            bc_f_mu_R_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(bc_f_mu_R_fusion_vertex_label))); #Hermitian Conjugate of the Splitting Vertex f -> bc through mu_R
                            af_d_nu_R_fusion_vertex_label := Concatenation(irrep_a, "*", irrep_f, "_", irrep_d, "__", String(nu_R));
                            af_d_nu_R_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(af_d_nu_R_fusion_vertex_label))); #Hermitian Conjugate of the Splitting Vertex d -> af through nu_R
                            #Compute diagram 2 via diagram stacking 
                            diagram2 := af_d_nu_R_fusion_vertex*KroneckerProduct(a_identity, bc_f_mu_R_fusion_vertex);
                            diagram2_norm_squared := [trivalent_vertices.(Concatenation("NORMSQ_", bc_f_mu_R_fusion_vertex_label)) , trivalent_vertices.(Concatenation("NORMSQ_", af_d_nu_R_fusion_vertex_label))];

                            #Add this left hand tree into the dictionary
                            rightTrees_record.(rightTree_label) := diagram2;
                            rightTrees_record.(rightTree_norm_squared_label) := diagram2_norm_squared;
                            #Add the key corresponding to this tree in the leftTrees_ordered_keys array
                            Add(rightTrees_ordered_keys, rightTree_label);

                        od;
                    od;
                od;

                #(3)-------------COMPUTATION OF GENERAL F SYMBOLS----------------------

                PrintTo(stream, "Here is the F matrix F^{", irrep_a, ",", irrep_b, ",", irrep_c, "}_", irrep_d, ": \n");
                
                for leftTree in leftTrees_ordered_keys do
                    for rightTree in rightTrees_ordered_keys do
                        #Create the F symbol label. "a,b,c,d,(e, mu_L, nu_L), (f, mu_R, nu_R)"
                        F_symbol_label := Concatenation(irrep_a, ",", irrep_b, ",", irrep_c, ",", irrep_d, ",(", leftTree, "),(", rightTree, ")");
                        
                        #Compute the F Symbol via Equation 2.8 (see diagram stacking version of Eq 2.8) 
                        diagram1 := rightTrees_record.(rightTree) ;
                        diagram2 := leftTrees_record.(leftTree) ;
                        diagram1_norm_squared := rightTrees_record.(Concatenation("NORMSQ_", rightTree));
                        diagram2_norm_squared := leftTrees_record.(Concatenation("NORMSQ_", leftTree));  

                        stacked_diagram := diagram2*TransposedMat(ComplexConjugate(diagram1));
                        #The F symbol is simply the trace of stacked_diagram divided by the dimension of irrep d
                        d_degree := DegreeOfRepresentation(irrep_dictionary.(irrep_d));
                        F_symbol_value := TraceMatrix(stacked_diagram)/d_degree;
                        
                        F_symbol_record.(F_symbol_label) := F_symbol_value;
                        F_symbol_norm_squared := diagram1_norm_squared[1]*diagram1_norm_squared[2]*diagram2_norm_squared[1]*diagram2_norm_squared[2];
                        

                        #OUTPUT THE DATA
                        PrintTo(stream, "(", leftTree, "),(", rightTree, ") = ");
                        PrintTo(stream, (F_symbol_record.(F_symbol_label)), "\n");
                        PrintTo(stream, "Normalization Squared Exact: ", F_symbol_norm_squared, "\n");
                        #PrintTo(stream, "Components: [", diagram1_norm_squared[1], ", ", diagram1_norm_squared[2], ", ", diagram2_norm_squared[1], ", ", diagram2_norm_squared[2], "] \n");
                        
                        PrintTo(stream, "\n \n");

                        #increment the number of non-trivial F symbols computed
                        total_non_trivial_F_symbols := total_non_trivial_F_symbols + 1;

                    od;
                od;
            od;
        od; 
    od;
od;

#Record GAP runtime when computation of F symbols ends
Fsymbol_computation_end := Runtime(); 

#Print Success Message & Report Runtime for Computation of F Symbols
Fsymbol_computation_time := Fsymbol_computation_end - Fsymbol_computation_start;
if Fsymbol_computation_time = 0 then    
    Print("\nAll F matrices for " , StructureDescription(G), " have successfully been computed! [< 1 ms] \n \n");
else   
    Print("\nAll F matrices for " , StructureDescription(G), " have successfully been computed! [", Fsymbol_computation_time," ms] \n \n");
fi;

CloseStream(stream);

#Report total runtime for the main algorithm
if Fsymbol_computation_end = intertwiner_computation_start then
    Print("[Main Algorithm Total Runtime: [< 1 ms] \n");
else
    Print("[Main Algorithm Total Runtime: ", Fsymbol_computation_end - intertwiner_computation_start ," ms] \n");
fi;
#===================================OUTPUT MESSAGES======================================
#Add in Runtime Data & F Symbol Overview into the Overview.txt file
filename:= Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Overview.txt");
stream := OutputTextFile(filename, true); #true - append mode on (since we are adding to the overview file here)

#Add F Symbols Overview
PrintTo(stream, "\n-----------Output Overview----------- \n");
PrintTo(stream, "Number of non-trivial F Matrices Computed: ", total_non_trivial_F_matrices ,"\n");
PrintTo(stream, "Number of non-trivial F Symbols Computed: ", total_non_trivial_F_symbols ,"\n \n");


#Add Runtime
PrintTo(stream, "-----------Runtime Breakdown--------- \n");
PrintTo(stream, "Note: A runtime of 0 ms indicates that the runtime was less than 1ms \n");
PrintTo(stream, "Time taken for Preprocessing: ", preprocessing_time_taken ," ms \n");
PrintTo(stream, "Time taken to compute full rank intertwiners: ", intertwiner_computation_time ," ms \n");
PrintTo(stream, "Time taken to orthogonalize columns of full rank intertwiners: ", orthogonalization_time_taken ," ms \n");
PrintTo(stream, "Time taken to extract trivalent vertices: ", trivalent_vertex_extraction_time ," ms \n");
PrintTo(stream, "Time taken to compute & store F symbols: ", Fsymbol_computation_time ," ms \n");
PrintTo(stream, "Total runtime: ", Fsymbol_computation_end - preprocessing_start  ," ms \n");

PrintTo(stream, "\n-----Statistics on Intertwiner Generation via Random Matrices----\n");
PrintTo(stream, "Entries of all random matrices were sampled uniformly from the range [-d,d]\\{0} where d = ", random_matrix_range, "\n");
PrintTo(stream, "Maximum number of re-runs required to produce a full rank intertwiner: ", max_random_matrix_runs[1], "\n");
if max_random_matrix_runs[1] > 1 then  
    PrintTo(stream, "The max number of runs were needed for the following ordered pairs of irreducible representations: ", max_random_matrix_runs[2], "\n");
fi;
CloseStream(stream);

Print("========================Output Data=================\n");
Print("The output has been stored in the following files: \n");
Print("(1) ", Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Overview.txt") ," - a general overview of Rep(", StructureDescription(G) ,") and its structure \n");
Print("(2) ", Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Output_Data.txt") ," - Contains all Fusion Rules and F matrices \n");
Print("(3) ", Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Unitary_Representations.txt") ," - Contains explicit unitary representations of ", StructureDescription(G) ," \n");
Print("======================END OF REPORT===========================\n \n");


#-----------------------------OUTPUT DATA IN FILES-----------------------------------
#-----------IGNORE THIS: NEEDS TO BE DEVELOPED FOR LATEX FILE OUTPUTS----------------- 

# filename := Concatenation("F_Symbols.txt", ".txt");
# stream := OutputTextFile(filename, false); #true - append mode on. false would mean create new file/empty existing one
# Print("The F symbol data has been printed to F_Symbols.txt ! \n");

# PrintTo(stream, "This file contains F Symbols for ", StructureDescription(G), "\n");



# for a in [2..repCount] do #pick first irrep label irrepNames[a]
#     for b in [2..repCount] do #pick second irrep label irrepNames[b]
#         for c in [2..repCount] do #pick third irrep label irrepNames[c] (here iterate over all non-trivial irreps)
#             #We are working with V^{abc} in this iteration of the loop
#             #We want to print F matrices for the fusion spaces V^{abc}_d for all valid d

#             #First we recall all valid d (i.e. access all the valid fusion outcomes of a,b and c via the triplet_fusion_record dictionary)
#             triplet_label := Concatenation(irrepNames[a], "*", irrepNames[b], "*", irrepNames[c]);
#             abc_distinct_outcomes := triplet_fusion_outcomes.(triplet_label); #access all the distinct fusion outcomes of a,b and c 
#             PrintTo(stream, triplet_fusion_outcomes, "\n");

#             #Now, we iterate over all distinct fusion outcomes of a,b and c V^{abc}_d for all d
#             for abc_outcome in abc_distinct_outcomes do
                
#                 #(3)-------------PRINTING GENERAL F SYMBOLS----------------------
#                 irrep_a :=  irrepNames[a];
#                 irrep_b := irrepNames[b];
#                 irrep_c := irrepNames[c];
#                 irrep_d := abc_outcome;
#                 PrintTo(stream, "Here is the F matrix F^{", irrep_a, ",", irrep_b, ",", irrep_c, "}_", irrep_d, ": \n");
#                 for leftTree in leftTrees_ordered_keys do
#                     for rightTree in rightTrees_ordered_keys do
#                         #Create the F symbol label. "a,b,c,d,(e, mu_L, nu_L), (f, mu_R, nu_R)"
                        
#                         F_symbol_label := Concatenation(irrep_a, ",", irrep_b, ",", irrep_c, ",", irrep_d, ",(", leftTree, "),(", rightTree, ")");
#                         PrintTo(stream, "(", leftTree, "),(", rightTree, ") = ");
#                         PrintTo(stream, F_symbol_record.(F_symbol_label));
#                         PrintTo(stream, "\n \n");

#                     od;
#                 od;
#             od;
#         od; 
#     od;
# od;
