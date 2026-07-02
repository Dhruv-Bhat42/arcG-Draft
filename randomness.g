#We define functions here to generate random integers, matrices and invertible matrices

#-------RANDOM INTEGER GENERATOR---------------
#We use GAP's in-built random number generator, the MersenneTwister. It always requires an initial seed value, which is hashed into a sequence of random numbers
#Choose the seed value to be internal runtime of GAP, so that the seed value is not hardcoded. 
#It will keep varying whenever this file is used. Otherwise, the same sequence of random numbers will be used over and over for different calls of this file
MyRandomSource := RandomSource(IsMersenneTwister, Runtime()); 

#Execute meaningless computations for a few microseconds to change the internal runtime of GAP
for futile_index1 in [1..10^6] do
    futile_index2 := futile_index1^2;
od;

# Define a function to generate a random integer between the bounds lower_bound and upper_bound (these must be integers, but NOT ZERO)
RandomInteger := function(lower_bound, upper_bound)
    local randNum;
    randNum := Random(MyRandomSource, lower_bound, upper_bound);
    #Make sure 0 is not produced
    while randNum = 0 do
        randNum := Random(MyRandomSource, lower_bound, upper_bound);
    od;
    return randNum; #is a random integer between -n and n
end;


#-------RANDOM MATRIX GENERATOR-----------------
#A function to generate a dim x dim matrix sampling from the set of integers in {lower_bound <= x <= upper_bound}
GenerateRandomMatrix := function(dim, lower_bound, upper_bound)
    local rand_matrix, row, col;
    rand_matrix := NullMat(dim, dim); #initializes a dim x dim zero matrix
    for row in [1..dim] do
        for col in [1..dim] do
            rand_matrix[row][col] := RandomInteger(lower_bound, upper_bound);
        od;
    od;
    return rand_matrix;
end;
