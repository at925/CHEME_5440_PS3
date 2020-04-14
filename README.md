## CHEME_5440_PS3

The Stoichiometric Matrix can be found in the ``Network.csv`` file, while the elemental composition of the metabolites that are being analyzed can be found in the ``Elements.csv`` file. 

**This information is presented in a more digestable form in the ``PS3.xlsx`` file detailing the reactions, stoichiometric matrix and elemental composition**

### Part A
The Stoichiometric Matrix consists for 18 rows (metabolites) and 21 columns (reactions and exchanges). The detailed breakdown can be found in ``PS3.xlsx`` while the matrix used in the Julia code can be found in ``Network.csv``

### Part B
To determine if the urea cycle reconstruction is elementally balanced for C,N,O,P and S (Elemental composition contained in ``Elements.csv``) issue the following command in Julia REPL:

 ```jl
    julia > include("Check_Element_Balance.jl")
  ```
  
This will compute the product of the transpose of the element matrix and the stoichiometric matrix resulting in a ``6x21 array``. We see that the first 6 columns (all the internal reactions) are 0 and are therefore balanced. 

### Part C
To estimate the maximum urea flux, issue the following command in Julia REPL:

```jl
    julia > include("Solve.jl")
  ```
Running the ``Solve`` script returns the objective_value which corresponds to the maximum urea flux. The ``calculated_flux_array``contains the optimal flux distributions for all metabolites in the case of optimal maximum urea production. The optimal urea flux calculated was found to be ``1.242 mmol/gDw-hr.``


  
