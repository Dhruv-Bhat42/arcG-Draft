#---------LOAD PACKAGES FOR CHAR TABLES, IRREPS AND UNITARISATION COMMANDS-------- 
if not IsPackageLoaded("ctbllib") then
    LoadPackage("ctbllib"); 
fi;

if not IsPackageLoaded("RepnDecomp") then 
    LoadPackage("RepnDecomp");   
fi;

#---------DEFINE GROUP, COMPUTE UNITARISED IRREDUCIBLE REPRESENTATIONS----------
#Computing Irreps in the sorted order of irreducible characters
Print("Computing Irreducible Representations of G via the Dixon Method... \n");
irreps := IrreducibleRepresentationsDixon(G, sortedChars);
Print("Success! \n");

#Unitarize Irreducible Representations
Print("Unitarizing Irreducible Representations for G... \n");
irreps_unitary := [];
iteration_num := 1;
for rho in irreps do
    if IsUnitaryRepresentation(rho) then
        Add(irreps_unitary, rho);
    else
        unitarized_rho:= UnitaryRepresentation(rho); 
        Add(irreps_unitary, unitarized_rho.unitary_rep);
        
    fi;
od;
Print("Unitarization Complete! \n \n");

# D10 := G;

# gens := GeneratorsOfGroup(G);

# if Order(gens[1]) = 5 then
#     sigma := gens[1];
#     tau := gens[2];
# else
#     sigma := gens[2];
#     tau := gens[1];
# fi;

# omega := E(5);
# zero := 0 * omega;

# piPlusSigma := [ [ 1 ] ];
# piPlusTau := [ [ 1 ] ];

# piPlus := GroupHomomorphismByImages(
# D10,
# Group(piPlusSigma, piPlusTau),
# [ sigma, tau ],
# [ piPlusSigma, piPlusTau ]
# );

# piMinusSigma := [ [ 1 ] ];
# piMinusTau := [ [ -1 ] ];

# piMinus := GroupHomomorphismByImages(
# D10,
# Group(piMinusSigma, piMinusTau),
# [ sigma, tau ],
# [ piMinusSigma, piMinusTau ]
# );

# pi1Sigma := [
# [ omega, zero ],
# [ zero, omega^-1 ]
# ];

# pi1Tau := [
# [ zero, 1 ],
# [ 1, zero ]
# ];

# pi1 := GroupHomomorphismByImages(
# D10,
# Group(pi1Sigma, pi1Tau),
# [ sigma, tau ],
# [ pi1Sigma, pi1Tau ]
# );


# pi2Sigma := [
# [ omega^2, zero ],
# [ zero, omega^-2 ]
# ];

# pi2Tau := [
# [ zero, 1 ],
# [ 1, zero ]
# ];

# pi2 := GroupHomomorphismByImages(
# D10,
# Group(pi2Sigma, pi2Tau),
# [ sigma, tau ],
# [ pi2Sigma, pi2Tau ]
# );

# irreps_unitary := [
# piPlus,
# piMinus,
# pi1,
# pi2
# ];



#------------------OUTPUT UNITARY IRREPS INTO A TEXT FILE-------------
# Open a file for writing the output
filename := Concatenation(String(IdGroup(G)[1]), ",", String(IdGroup(G)[2]), "_Unitary_Representations.txt");
stream := OutputTextFile(filename, false); #true - append mode/ false - create new file each time
PrintTo(stream, "This file contains explicit (unitary) irreducible representations for ", StructureDescription(G),  " \n");
PrintTo(stream, "The GAP I.D. for this group is ", IdGroup(G), "\n");

# Write unitarized irreducible representations into the text file
for i in [1..Length(irreps_unitary)] do
    PrintTo(stream, "Irreducible Representation ", i, ":\n");
    
    
    # Apply the representation to each group element
    for g in Elements(G) do
    	PrintTo(stream, "Group Element: ", g, ": \n");
        PrintTo(stream, irreps[i](g), "\n");

    od;
    
    PrintTo(stream, "\n");  # Add a newline for better separation
od;

# Close the file stream
CloseStream(stream);
