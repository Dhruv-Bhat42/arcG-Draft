
#Designing a subroutine to separate multiplicities from irrep labels. 
#For example input may be "2R_3". So we need to extract multiplicity 2 and the irrep "R_3"
label_info := function(irrep_label) 
#Declare local variables
local multiplicities_present, multiplicity_string_index, current_multiplicity, current_irrep_label, label_info;
irrep_label := String(irrep_label);

#First we check if multiplicities are contained
multiplicities_present := false;
if Position(irrep_label, 'R') > 1 then
    multiplicities_present := true;
fi;

#Now we will create and return a list [n, R_d_j] where n is the multiplicity of "irrep_label" and R_d_j is the irrep label contained within "irrep_label"
if multiplicities_present then
    
    multiplicity_string_index := Position(irrep_label, 'R') - 1; #This returns the [(index of R) - 1] in the irrep label string. The multiplicity is stored from position 1 to this position.                    
    #The multiplicity of the current irrep is the irrep_label string, sliced from position 1 to position mulitplicity_string_index
    current_multiplicity := Int(irrep_label{[1..multiplicity_string_index]});
    #The irrep label itself that irrep_label contains
    current_irrep_label := irrep_label{[(multiplicity_string_index+1)..Length(irrep_label)]};

    label_info := [current_multiplicity, current_irrep_label];

else
    label_info := [1, irrep_label];
fi;
   
return label_info;
end;

