# arc(G) User Guide

## Description

“arc(G)” is a software tool written in GAP that can be used to compute the F symbols (associators) for representation categories of finite groups. This is an end-user guide for individuals attempting to use arcG.

Via arcG, you have the ability to:

- Select any finite group of your choosing (of reasonably low order)
- View the fusion ring associated to this group, and learn about key properties such as the presence of multiplicities, the number & degrees of irreducible representations.
- Gain access to irreducible unitary representations of this group (this uses the existing GAP package “Repspn”)
- Gain access to F symbols of the representation category of the group in unitary gauge, and certain statistics related to them (total number of non-trivial F matrices and non-trivial F symbols)
- Modify certain algorithm parameters to make a unitary change of gauge for simplifications

## Section 1: GAP installation and Setup for arcG

### Part A (Installations)

1. You must first install GAP on your laptop/PC. Visit <https://www.gap-system.org/install/> and follow the particular instructions for your Windows/MAC/Linux system.

2. (Optional) Use the instructions on the webpage to test that your installation was successful.

3. DVisit our GitHub and download “arcG.zip” your device. Unzip/Extract this zip file into a folder at a suitable location and name it “arcG”.

### Part B (Setting up arcG)

1. **Launch the GAP terminal and note the current working directory:** Launch the GAP terminal, type in

    ```gap
    DirectoryCurrent();
    ```

    and press enter. Currently, this is the folder from which the terminal can read “.g” files (GAP scripts).

    Example:

    <!-- Insert screenshot/example image here if desired. -->

2. **Changing the working directory to the arcG folder:** Copy the file path of the arcG folder and in the GAP terminal, type

    ```gap
    ChangeDirectoryCurrent();
    ```

    and paste the local path of the arcG folder inside the parenthesis in the same syntax as you observe in the output of

    ```gap
    DirectoryCurrent();
    ```

    For example, since my arcG folder is saved at Desktop, then, the correct command would be

    ```gap
    ChangeDirectoryCurrent("/cygdrive/c/Users/Dhruv/Desktop/arcG");
    ```

    **Warning:** Note that the correct GAP syntax for the path may differ from how you see it on your device. For example, as a Windows user, I see the file path as in

    <!-- Insert Windows filepath screenshot/example image here if desired. -->

    But the GAP syntax starts with `cygdrive/c` and uses forward slashes instead of backward slashes.

3. **Test that you are in the appropriate working directory:** If the terminal doesn’t return any errors in Step 4, type

    ```gap
    DirectoryCurrent();
    ```

    and press enter. This should now reflect the path for the arcG folder. If this works, the setup is complete!

4. **Create a change-of-directory shortcut:** GAP will revert back to the original directory (as seen in Step 1) every time you close the terminal and relaunch it. For this purpose, it’s best to save the step 2 command somewhere. Here is a simple example Open any text-editing system (like notepad) and type in the correct command from step 2. Save this file as “ChangePath.g” in the folder that the terminal always resets to (from Step 1). Now, every time you launch the GAP terminal, you need only type in

    ```gap
    Read("ChangePath.g");
    ```

## Section 2: Using arcG

After you have navigated to the appropriate directory corresponding to arcG (see section 1 part B), you are ready to interact with arcG.

Type

```gap
Read("main.g");
```

and hit enter to begin interacting with arcG.

### Step 1

You will immediately be prompted to input the name of the finite group you want to work with. You can learn about the GAP naming conventions here. GAP assigns an ID of form $[n,j]$ to all finite groups of sufficiently low order, wherein $n$ represents the order of the group and $j$ is an arbitrary indexing within isomorphism classes of groups of that order. Here are some simple examples of commands based on GAP’s library of finite groups:

| Command | Explanation |
|---|---|
| `G := SymmetricGroup(3);` | Assigns the group $S_3$ to the variable `G` |
| `G := AlternatingGroup(4);` | Assigns the group $A_4$ to the variable `G` |
| `G := DihedralGroup(n);` | Assigns the Dihedral group of order $n$ to the variable `G`. Since $n$ must always be even, returns an error if $n$ is odd. |
| `G := SmallGroup(n,j);` | For a particular choice of positive integers $n$ and $j$, the command assigns a particular finite group of order $n$ to the variable `G` (as per GAP’s internal indexing). If $j$ is greater than the total number of isomorphism classes of order-$n$ groups, then the command returns an error. |

After you assign a finite group to `G`, you can simply press enter.

### Step 2

You will be asked to input a positive integer $d$ between 1 and 10. At a certain point in the algorithm, we use random integer matrices wherein entries are sampled uniformly from the set $[-d,d]\setminus\{0\}$. If $d$ is small, the probability that a specific step in the algorithm- which relies on random matrices- fails and needs to be rerun is greater. But the upside is that the final F symbols are computed in a simpler gauge. If $d$ is higher value, the script may take very long to run and the gauge for the F symbols is decidedly less visually palatable. We highly advise that you choose $d = 1$ for simplicity and speed.

When working with Abelian groups, $d=1$ is always the best choice, since the step involving random matrices never fails in this case.

### Potential Failure

After Step 2, the script will run until it produces all the final output, or runs into an error. One of the main reasons for an error may be that GAP is unable to unitarize irreducible representations for the user’s choice of group. The unitarization for the irreducible representations of these groups involves handling nested radicals, for which GAP is unable to process exact computations. Section 4 contains a list of finite groups of order less than or equal to 32 for which GAP is unable to unitarize representations using the in-built unitarization command.

### Output

If there are no errors, the output will consist of the following files. Note that “n,j” represents the GAP ID of the input group

1. An Overview File titled `n,j_Overview.txt` containing the name of the group `G`, the GAP ID of `G`, total number of irreducible representations, whether multiplicities are present, irreducible representation labels, fusion rules, and the number of non-trivial F matrices & F symbols computed. A list of runtimes for various phases of the script such as Preprocessing, Computing Full Rank Intertwiners, Orthogonalization, F Symbols, and Total Compute Time are presented.

2. List of Unitarized Irreducible Representations of `G` in a file titled `n,j_Unitary_Representations.txt`

3. Main Output File containing all F symbol data of `G` (with informative headers) in a file titled `n,j_Output_Data.txt`. See section 3 for its interpretation.

## Section 3: Data Interpretation

### Irreducible Representation Labels

An isomorphism class of irreducible representations for a finite group is labelled as $R_d$ where $d$ represents the dimension of the representation. If there are more than one isomorphism classes of irreducible representations of dimension $d$, a further (arbitrary) indexing $R_{d_j}$ will be used to distinguish between them. $R_I$ is used to refer to the trivial irreducible representation.

### Sample Overview File Example

The alternating group $A_4$ consists of two non-trivial isomorphism classes of irreducible representations of degree 1 and one irreducible representation of degree 3. Thus, it has the label set $\{R_I,R_{1_2},R_{1_3},R_3\}$. Moreover, $A_4$ has GAP ID $[12,3]$. An arcG output file for $A_4$ would resemble this text file:

<!-- Insert sample overview file screenshot or text block here if desired. -->

Note that $\operatorname{Rep}(A_4)$ is not multiplicity-free, which is reported in the 4th line (and also in the last non-trivial fusion rule $R_3*R_3=R_I+R_{1_2}+R_{1_3}+2R_3$)

### Main Output File Interpretation

Each F matrix is characterized by 4 labels. arcG sequentially prints non-trivial F matrices under the headline

```text
Here are the F symbols for the matrix F_d^abc
```

The F matrices are printed in ascending order of representation degree in each index of the F matrix from left to right. For example, for $\operatorname{Rep}(A_4)$, the first matrix printed is $F_{R_I}^{R_{1_2},R_{1_2},R_{1_2}}$. Since $R_I$ is the only fusion outcome for this label triplet, the variable c will be incremented next, and the next matrix would be $F_{R_{1_2}}^{R_{1_2},R_{1_2},R_{1_3}}$ ($R_{1_2}$ is the unique fusion outcome for this triplet). We continue this process until c and d take on all possible values in ascending order, and then increment the label b, and the process continues.

### F symbol labels

The F symbol $[F_d^{abc}]_{(e_{\mu_\nu},f_{\mu'_{\nu'}})}$ is labelled simply as the tuple $(e_{\mu_\nu},f_{\mu'_{\nu'}})$ under the $F_d^{abc}$ heading. These are also printed in ascending order of representation degree. [Can add sample matrix created by Sachin here]

For each F symbol, the data is printed in two parts:

a. F symbol numerator (computed via diagram stacking)

b. F symbol denominator squared (computed by multiplying normalizations of every trivalent vertex involved in diagram stacking)

Thus, the F symbol $[F_d^{abc}]_{(e_{\mu_\nu},f_{\mu'_{\nu'}})}$ should be interpreted as $\frac{a}{\sqrt{b}}$ where a and b are as above.

In case the representation category contains no multiplicities, the F symbol labels still contain the indices $\mu,\nu,\mu'$ & $\nu'$, but each of these are set to 1, so they may be ignored.

### Sample Output  \(\operatorname{Rep}(A_4)\)

[Add pictures of GAP output and LaTeX output for comparison]

## Section 4: Miscellaneous Details

### Failure of Unitarization for Irreducible Representations

Here is a list of finite groups of order less than or equal to 32 for which GAP is unable to unitarize representations using the in-built unitarization command.

An entry $[n,j]$ in the following array corresponds to the `GAP SmallGroup(n,j)`.

```text
[[10,1], [14,1], [18, 1], [20,1], [20,4], [22,1], [26,1], [28,1], [28,3], [30,2], [30,3], [32,18], [32,19], [32,20]]
```

### Statistics on Intertwiner Generation via Random Matrices

For range $[-d,d]\setminus\{0\}$

“Over X trials, the max number of re-runs of the algorithm needed for a full rank intertwiner was N, for all finite groups of order < Y”

| $d$ (Integer Input Range) | X (number of trials) | N (max re-runs) | Y (over all groups of order) |
|---:|---:|---:|---:|
| 1 | 30 | 3 | 38 |
| 5 | 20 | 2 | 32 |
| 10 | 20 | 2 | 32 |

### Tips

Don’t forget to add the semicolon after each command

Error-handling: If you come across an error, type

```gap
quit;
```

and hit enter to escape the break loop.

GAP indexing starts at 1, not 0.
