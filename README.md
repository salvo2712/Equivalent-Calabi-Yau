# Equivalent-Calabi-Yau
In this code, we will use CYTools to see if two CY from different polytopes, but with the same Hodge numbers are topologically equivalent or not, following the Wall's theorem.

## Description

For each polytope retrieved via CYTools, the code:
- Computes a triangulation and constructs the associated Calabi-Yau threefold
- Extracts the Hodge numbers (h11, h21), the triple intersection numbers κijk, 
  and the second Chern class c2
- Groups CY manifolds by Hodge pair (h11, h21)
- For each pair of CY manifolds sharing the same Hodge numbers, searches for an 
  integer basis change matrix M ∈ GL(h11, ℤ) such that:
  - M · c2(A) = c2(B)
  - The intersection tensor transforms as κ'_ijk = Σ M_ix M_jy M_kz κ_xyz

If such a matrix exists, the two manifolds are identified as equivalent 
(i.e. related by a change of basis in H2).

## Requirements

- Python 3
- [CYTools](https://cy.tools)
- NumPy

## Usage

Run the script directly or open it in JupyterLab:
```bash
python EQUIVALENT_CY.py
```

The default analysis targets polytopes with `h11=2`. The search range for 
matrix entries can be tuned via the `limit` parameter in 
`trova_matrice_equivalenza()`.

## Output

The script prints all pairs of CY manifolds found to be equivalent, along with 
the explicit basis change matrix and a verification of the c2 transformation.
