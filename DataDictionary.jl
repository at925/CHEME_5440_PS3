using CSV
using DataFrames
using LinearAlgebra

function build_data_dictionary()

    data_dictionary = Dict{AbstractString,Any}()

    stoichiometric_matrix = DataFrame(CSV.File("Network.csv",header=false))
    stoichiometric_matrix = convert(Matrix,stoichiometric_matrix)
    data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix

    # Enter specific values for v1 -> v5 reactions

    E = 0.01E-3 # steady-state enzyme concentration (units [mmol/gDW])

    # Metabolic kcats for v1 -> v5
    kcat_v1 = (203)*(3600)   # v1 ATP + L-Citrulline + L-Aspartate --> AMP + Diphosphate + N-(L-Arginino)succinate (units [1/hr])
    kcat_v2 = (34.5)*(3600)  # v2 N-(L-Arginino)succinate --> Fumarate + L-Arginine (units [1/hr])
    kcat_v3 = (249)*(3600)   # v3 L-Arginine + H2O --> L-Ornithine + Urea (units [1/hr])
    kcat_v4 = (88.1)*(3600)  # v4 Carbamoyl_phosphate + L-Ornithine --> Orthophosphate + L-Citrulline (units [1/hr])
    kcat_v5f = (13.7)*(3600) # v5f 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H --> 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O (units [1/hr])
    kcat_v5r = (13.7)*(3600) # v5r 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O --> 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H (units [1/hr])

    # Metabolic V-max for v1 -> v5
    v1_max = (kcat_v1)*(E)
    v2_max = (kcat_v2)*(E)
    v3_max = (kcat_v3)*(E)
    v4_max = (kcat_v4)*(E)
    v5f_max = (kcat_v5f)*(E)
    v5r_max = (kcat_v5r)*(E)

    # Cell Dry Weight Calculation
    cell_mass = (2.3E-9)    # Estimated HeLa cell mass from BioNumbers (ID 103720) (units g)
    cell_volume = 2.425E-12 # HeLa cell volume from BioNumbers (ID 103725) (units L)
    water_fraction = 0.70   # Estimated water ratio in HeLa cells from BioNumbers (ID 100387)
    cell_dry_mass = (1-water_fraction)*cell_mass
    cell_dry_weight = (cell_volume/cell_dry_mass)

    #Convert Metabolite concentration and Km values:-
    ATP_conc = (4.67E-3*1000)*(cell_dry_weight)       # v1 for Mus musculus from Park et. al
    ATP_Km = (3.92E-4*1000)*(cell_dry_weight)         # v1 for Mus musculus from Park et. al
    aspartate_conc = (1.49E-2*1000)*(cell_dry_weight) # v1 for Mus musculus from Park et. al
    aspartate_Km = (1.54E-4*1000)*(cell_dry_weight)   # v1 for Mus musculus from Park et. al
    arginine_conc = (2.55E-4*1000)*(cell_dry_weight)  # v3 for Homo sapiens from Park et. al
    arginine_Km_1 = (1.55E-3*1000)*(cell_dry_weight)  # v3 for Homo sapiens from Park et. al
    arginine_Km_2 = (3.50E-6*1000)*(cell_dry_weight)  # v5f for Mus musculus from Park et. al
    ornithine_conc = (4.49E-3*1000)*(cell_dry_weight) # v4 for Saccharomyces cerevisiae from Park et. al
    ornithine_Km = (1.60E-3*1000)*(cell_dry_weight)   # v4 for Saccharomyces cerevisiae from Park et. al

    # Metabolic V for v1 -> v5
    v_1 = v1_max*((ATP_conc)/(ATP_Km + ATP_conc))*((aspartate_conc)/(aspartate_Km + aspartate_conc))
    v_2 = v2_max*(1) # As no metabolite concentration/Km could be found for human cell line, saturation term = 1
    v_3 = v3_max*((arginine_conc)/(arginine_Km_1 + arginine_conc))
    v_4 = v4_max*(1) # As no metabolite concentration/Km could be found for human cell line, saturation term = 1 (Note the ornithin numbers fonud were for Saccharomyces cerevisiae
    v_5f = v5f_max*((arginine_conc)/(arginine_Km_2 + arginine_conc))
    v_5r = v5r_max*(1) # As no metabolite concentration/Km could be found for human cell line, saturation term = 1

    # Building the Metabolic Flux Bounds Array
    metabolic_flux_bounds_array = [
    0.0 v_1;       # v1 (units [mmol/gDW-hr])
    0.0 v_2;       # v2 (units [mmol/gDW-hr])
    0.0 v_3;       # v3 (units [mmol/gDW-hr])
    0.0 v_4;       # v4 (units [mmol/gDW-hr])
    0.0 v_5f;      # v5f (units [mmol/gDW-hr])
    -v_5r 0.0 ;    # v5r (units [mmol/gDW-hr])
    0.0 10.0;           # b1 [] -> Carbamoyl_phosphate (units [mmol/gDW-hr])
    0.0 10.0;           # b2 [] -> L-Aspartate (units [mmol/gDW-hr])
    0.0 10.0;           # b3 Fumarate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b4 Urea -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b5 [] -> ATP (units [mmol/gDW-hr])
    0.0 10.0;           # b6 AMP -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b7 Diphosphate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b8 Orthophosphate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b9 [] -> Oxygen (units [mmol/gDW-hr])
    0.0 10.0;           # b10 [] -> NADPH (units [mmol/gDW-hr])
    0.0 10.0;           # b11 [] -> H (units [mmol/gDW-hr])
    0.0 10.0;           # b12 Nitric_oxide -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b13 NADP -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b14 [] -> H20 (units [mmol/gDW-hr])
    0.0 10.0;           # b15 H20 -> [] (units [mmol/gDW-hr])
    ]

    # Setup default flux bounds
    data_dictionary["metabolic_flux_bounds_array"] = metabolic_flux_bounds_array
    # Setup default species bounds array
    species_bounds_array = [
    0.0 0.0; #1 ATP
    0.0 0.0; #2 L-Citrulline
    0.0 0.0; #3 L-Aspartate
    0.0 0.0; #4 AMP
    0.0 0.0; #5 Diphosphate
    0.0 0.0; #6 N-(L-Arginino)succinate
    0.0 0.0; #7 Fumarate
    0.0 0.0; #8 L-Arginine
    0.0 0.0; #9 H20
    0.0 0.0; #10 L-Ornithine
    0.0 0.0; #11 Urea
    0.0 0.0; #12 Carbamoyl_phosphate
    0.0 0.0; #13 Orthophosphate
    0.0 0.0; #14 Oxygen
    0.0 0.0; #15 NADPH
    0.0 0.0; #16 H
    0.0 0.0; #17 Nitric_oxide
    0.0 0.0; #18 NADP
    ]

    data_dictionary["species_bounds_array"] = species_bounds_array

    # Setup the objective coefficient array
    objective_coefficient_array = [
    0.0;       # v1 (units [mmol/gDW-hr])
    0.0;       # v2 (units [mmol/gDW-hr])
    0.0;       # v3 (units [mmol/gDW-hr])
    0.0;       # v4 (units [mmol/gDW-hr])
    0.0;       # v5f (units [mmol/gDW-hr])
    0.0;       # v5r (units [mmol/gDW-hr])
    0.0;           # b1 [] -> Carbamoyl_phosphate (units [mmol/gDW-hr])
    0.0;           # b2 [] -> L-Aspartate (units [mmol/gDW-hr])
    0.0;           # b3 Fumarate -> [] (units [mmol/gDW-hr])
    1.0;           # b4 Urea -> [] (units [mmol/gDW-hr])
    0.0;           # b5 [] -> ATP (units [mmol/gDW-hr])
    0.0;           # b6 AMP -> [] (units [mmol/gDW-hr])
    0.0;           # b7 Diphosphate -> [] (units [mmol/gDW-hr])
    0.0;           # b8 Orthophosphate -> [] (units [mmol/gDW-hr])
    0.0;           # b9 [] -> Oxygen (units [mmol/gDW-hr])
    0.0;           # b10 [] -> NADPH (units [mmol/gDW-hr])
    0.0;           # b11 [] -> H (units [mmol/gDW-hr])
    0.0;           # b12 Nitric_oxide -> [] (units [mmol/gDW-hr])
    0.0;           # b13 NADP -> [] (units [mmol/gDW-hr])
    0.0;           # b14 [] -> H20 (units [mmol/gDW-hr])
    0.0;           # b15 H20 -> [] (units [mmol/gDW-hr])
    ]

    data_dictionary["objective_coefficient_array"] = objective_coefficient_array

    data_dictionary["min_flag"] = false

    return data_dictionary

end
