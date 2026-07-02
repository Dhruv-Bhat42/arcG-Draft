# `G_fusion_finder.g`
This page documents the GAP file `G_fusion_finder.g`.

## Mathematical Overview

The purpose of `G_fusion_finder.g` is to compute the nontrivial fusion rules of the representation category \(\operatorname{Rep}(G)\), where \(G\) is a finite group. 

For irreducible representations \(V_a\), \(V_b\) and \(V_c\), with corresponding characters \(\chi_a\), \(\chi_b\) and \(\chi_c\)  the fusion coefficient \(N_{c}^{ab} \) is defined by

$$
N_{c}^{ab} = \dim (\operatorname{Hom}_{Rep(G)}(V_a \otimes V_b, V_c)) = \langle \chi_a \cdot \chi_b , \chi_c \rangle
\label{eq:fusion-coefficient}
\tag{1}
$$

where \( \langle \cdot , \cdot \rangle \) is the inner product on class functions on G, defined by 

$$
\langle \chi,\psi\rangle
=
\frac{1}{|G|}
\sum_{g\in G}
\chi(g)\overline{\psi(g)}.
$$

The file `G_fusion_finder.g` computes these coefficients using the character table of \(G\), stores them into a variable called `fusionRules` and prints them in a text file (called `n_j_Overview.txt` where `(n,j)` is the GAP ID of \(G\))

---

## File Dependencies

- The file expects a finite group variable `G` to be defined before the main computation is carried out. This definition is inherited from `main.g` when the user is prompted to define a group.

- The file loads the character table of \(G\) using functions from GAP's `ctbllib` package. 

---

## Overview of the Script

`G_fusion_finder.g` performs the following steps:

1. Loads the `ctbllib` package.
2. Uses the `ctbllib` package to access the explicit irreducible characters of \(G\).
3. Sorts the irreducible characters of \(G\) by degree.
4. Assigns readable labels to each irreducible representation.
5. Computes fusion rules (the \( N_c^{ab} \) values)  using the character inner product. 
6. Stores the resulting fusion rules in a variable called `fusionRules`. 
7. Prints the nontrivial fusion rules in the text file (called `n_j_Overview.txt` where `(n,j)` is the GAP ID of \(G\))

---

## Extracting Irreducible Characters

The file extracts the irreducilbe characters of \(G\) using in-built commands from the `ctbllib` package. 

```gap
charTable := CharacterTable(G);
irreducibleCharacters := Irr(charTable);
```

The variable `charTable` stores the character table of the group \(G\).

The variable `irreducibleCharacters` stores the list of irreducible characters associated to this character table.

The file then performs a simple validation on the number of entries in `irreducibleCharacters`  based on the standard fact that the number of irreducible characters of \(G\) agrees with the number of conjugacy classes of \(G\):

```gap
repCount := Length(ConjugacyClasses(G));

if not repCount = Length(irreducibleCharacters) then
    Print("WARNING: Character Table was not loaded correctly \n");
fi;
```

---

## Sorting irreducible characters

The irreducible characters are sorted by their degree (the dimension of the corresponding irreducible representation) using the `SortedCharacters()` command. Sorting by degree makes it easier to label irreducible representations and make fusion rules easier to read. 

```gap
sortedChars := SortedCharacters(charTable, irreducibleCharacters, "degree");
```

The trivial representation always appear first in this sorting. If two irreducible characters have the same degree, this command arbitrarily assigns an order between them. This order is not altered for the remainder of the file. 

The file then records the degree of each irreducible character in the variable `irrepDegree`

```gap
irrepDegree := [];

for i in [1..Length(sortedChars)] do
    degree := Degree(sortedChars[i]);
    Add(irrepDegree, degree);
od;
```

The list `irrepDegree` stores the dimensions of the irreducible representations in the same order as `sortedChars`. Thus, the integer accessed by `irrepDegree[i]` is the dimension of the irreducible representation whose character is `sortedChars[i]`.

---

## Naming irreducible representations

Next, the file assigns readable labels to each irreducible representations.

The trivial representation is labeled \(R_{\mathrm{I}}\) in the program output. 

A nontrivial irreducible representation of dimension \(d\) is labeled \(R_d\).

If there is more than one irreducible representation of dimension \(d\), a secondary index is added as per the order in `sortedChars`. These representations are labeled

$$
R_{d,1}, R_{d,2}, R_{d,3}, \ldots
$$

For example, if a group has two irreducible representations of dimension \(3\), they are labeled

$$
R_{3,1}, R_{3,2}.
$$

In the GAP output, these same labels appear as `R_3_1` and `R_3_2`.

The labels are stored in the list

```gap
irrepNames := [];
```

The list `irrepNames` is ordered consistently with `sortedChars`. Therefore, `irrepNames[i]` is the label attached to the irreducible character `sortedChars[i]`.

The relevant section of the code is:

```gap
countWithinDimension := 1;
for i in [1..Length(sortedChars)] do
    
    #Labelling Trivial Irrep
    if i=1 then 
    	Add(irrepNames, "R_I"); #label of the trivial irrep
    	
   	else
    	irrepLabel := Concatenation("R_", String(irrepDegree[i])); #Label the irrep by the dimension
    	
    	#Handle the indexing for irreps of the same degree
    	if Degree(sortedChars[i-1]) < irrepDegree[i] then
    		countWithinDimension := 1; #first irrep of this dimension found
    	else
    		countWithinDimension := countWithinDimension + 1; #increment the number of irreps this dimension
    	fi;

    	#Add further indexing to irrep label if multiple irreps of the same degree exist
    	if countWithinDimension > 1 then
    		irrepLabel := Concatenation(irrepLabel, "_", String(countWithinDimension));
    	else
    		#If countWithinDimension is 1, but there are more to follow
    		if IsBound(sortedChars[i+1]) and Degree(sortedChars[i+1]) = irrepDegree[i] then
    			irrepLabel := Concatenation(irrepLabel, "_", String(countWithinDimension));
    		fi;
    	fi;

    	Add(irrepNames, irrepLabel); #Add this label to the list
    fi;
od;

```

---

## Fusion rule computation

The fusion rule computation is performed by looping over pairs of nontrivial irreducible characters in the following code block:

```gap
fusionRules := [];
multiplicities_present := false;

for i in [2..Length(sortedChars)] do
    for j in [i..Length(sortedChars)] do
        chi := sortedChars[i];
        psi := sortedChars[j];

        product := chi * psi;

        decomposition := [];

        for k in [1..Length(sortedChars)] do
            phi := sortedChars[k];
            coefficient := ScalarProduct(charTable, product, phi);

            if coefficient > 0 then
                if coefficient = 1 then
                    Add(decomposition, irrepNames[k]);
                else
                    multiplicities_present := true;
                    Add(decomposition, Concatenation(String(coefficient), irrepNames[k]));
                fi;
            fi;
        od;

        Add(fusionRules, decomposition);
    od;
od;
```


### Loop Indices
The outer loops start at index `2` because the first irreducible character corresponds to the trivial representation. The trivial representation acts as the monoidal unit, so tensor products involving it don't need to be computed:

\[
R_I \otimes R_d \cong R_d.
\]

The second loop starts at `j = i` because tensor products in \(\operatorname{Rep}(G)\) are commutative up to isomorphism. Thus, the file computes fusion rules only once per representative of each unordered pair.

### Fusion Rule Computation
For each pair of irreducible characters \(\chi\) and \(\psi\), the product

```gap
product := chi * psi;
```

represents the character of the tensor product representation. 

The code then finds each fusion coefficient by computing the inner product with each irreducible character, as in Equation \(\eqref{eq:fusion-coefficient}\):

```gap
coefficient := ScalarProduct(charTable, product, phi);
```

If the coefficient is \(1\), the code stores only the representation label:

```gap
Add(decomposition, irrepNames[k]);
```

If the coefficient is greater than \(1\), the code stores the coefficient together with the representation label, and flips the `multiplicities_present` variable to true:

```gap
Add(decomposition, Concatenation(String(coefficient), irrepNames[k]));
```

### Storing Fusion Rules
The fusion rules are stored in the variable `fusionRules`, which is used very frequently in ``main.g``.

```gap
fusionRules := []; #a list of lists
```

Each entry of `fusionRules` is a list describing the decomposition of one tensor product.

For example, if

\[
R_a \otimes R_b \cong R_I \oplus R_2 \oplus 2R_3.
\]

then this information is stored in ``fusionRules`` as the list

```gap
[ "R_I", "R_2", "2R_3" ]
```
The order of entries in `fusionRules` follows the nested loop

```gap
for i in [2..Length(sortedChars)] do
    for j in [i..Length(sortedChars)] do
        ...
    od;
od;
```

Therefore, the list `fusionRules` stores decompositions in an ascending order by irrep degree, avoiding any duplicates and trivial fusion rules. 

The advantage of this data structure is that it allows us to access all fusion rules sequentially, which we need to do multiple times later on. However, the disadvantage is that it is hard to *look up* the fusion rule for a particular pair. We eventually also need to do this, so we later create the dictionary data structure ``fusionRules_record`` in ``main.g`` for this purpose. 

---

## Outputting Data

This last section of ``G_fusion_finder.g`` creates two text files containing information about \(\operatorname{Rep}(G)\): 
1. `n,j_Output_Data.txt`, which is the main output-data file (into which F symbols are later printed); 
2. `n,j_Overview.txt`, which contains a summary of the group, its irreducible representations, and its nontrivial fusion rules. 

Here, \([n,j]\) is the GAP ID of the input group \(G\).

### Output to the main n,j_Output_Data.txt file 

The first block creates the file `n,j_Output_Data.txt` for the first time: 

```gap 
filename:= Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Output_Data.txt"); 
stream := OutputTextFile(filename, false); 
``` 

The script then writes a header containing a readable structural description of the group: 
```gap 
PrintTo(stream, "This file contains output data for Rep(G) where G = ", StructureDescription(G), "\n"); 
``` 

The GAP command ``StructureDescription(G) `` returns a human-readable description of the group \(G\), when GAP is able to provide one. For example, this may print a header resembling: 
```text 
This file contains output data for Rep(G) where G = A4 
``` 

The code then prints the list of irreducible-representation labels: 

```gap 
PrintTo(stream, "Irrep Labels: "); 
for i in [1..Length(irrepNames)] do 
    if i = Length(irrepNames) then 
        PrintTo(stream, irrepNames[i], "\n"); 
    else 
        PrintTo(stream, irrepNames[i], ", "); 
    fi; 
od; 
```

The list `irrepNames` was constructed earlier in the program and is ordered consistently with `sortedChars`. The loop prints each label followed by a comma, except for the final label, which is followed by a newline. Thus, the output has the form: 
```text 
Irrep Labels: R_I, R_1_1, R_1_2, R_3 
```

After writing this initial information, the stream is closed: 
```gap 
CloseStream(stream); 
``` 
Closing the stream is necessary before the file is later reopened to write the remaining output data.

### Output to n,j_Overview.txt

The last section of the file prints the nontrivial fusion rules into the text file ``n,j_Overview.txt`` where ``[n,j]`` is the GAP ID of the group \(G\). 

The next block creates the overview file: 
```gap 
filename:= Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Overview.txt"); stream := OutputTextFile(filename, false); 
```
This text file begins with two identifying lines: 
```gap 
PrintTo(stream, "This file contains an overview of Rep(G) where G = ", StructureDescription(G), "\n"); 
PrintTo(stream, "The GAP I.D. for this group is ", IdGroup(G), "\n"); 
```
And then displays two important properties of \(\operatorname{Rep}(G)\)- the rank and the presence of multiplicities using the variables ``repCount`` and ``multiplicities_present``: 
```gap 
PrintTo(stream, "The rank of Rep(G) (i.e. the number of isomorphism classes of irreducible representations) is ", repCount, "\n"); 
PrintTo(stream, "Multiplicity-free? : ", not multiplicities_present, "\n \n"); 
```
The overview file also displays all the irrep labels, just as in ``n,j_Output_Data.txt``. 

Finally, the following loop writes the previously computed entries of fusionRules to the overview file.

```gap
PrintTo(stream, "\n------Non-trivial Fusion Rules------\n");
repCount := Length(sortedChars);
fusionIndex := 1;

for i in [2..repCount] do
    for j in [i..repCount] do
        currentDecomp := fusionRules[fusionIndex];

        PrintTo(stream, irrepNames[i], " * ", irrepNames[j], " = ");
        for k in [1..Length(currentDecomp)] do
            if k = Length(currentDecomp) then
                PrintTo(stream, currentDecomp[k], "\n");
            else
                PrintTo(stream, currentDecomp[k], " + ");
            fi;
        od;

        fusionIndex := fusionIndex + 1;
    od;
od;
```

Because ``fusionRules`` was populated in the same nested-loop order, ``fusionIndex`` identifies the decomposition corresponding to the current pair ``irrepNames[i]`` and ``irrepNames[j]``.

The inner loop formats each stored decomposition as a sum of labels, producing output such as:

```text
R_a * R_b = R_c + R_d + ... 
```

