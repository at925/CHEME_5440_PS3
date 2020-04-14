using CSV
using LinearAlgebra
using DataFrames

stoichiometric_matrix = DataFrame(CSV.File("Network.csv",header=false))
stoichiometric_matrix = convert(Matrix,stoichiometric_matrix)
element_matrix = DataFrame(CSV.File("Elements.csv",header=false))
element_matrix = convert(Matrix,element_matrix)

element_balance = transpose(element_matrix )*stoichiometric_matrix
