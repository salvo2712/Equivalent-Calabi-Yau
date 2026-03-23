import numpy as np
from itertools import product, combinations
from cytools import fetch_polytopes


# ─────────────────────────────────────────────
# CELLA 1 — Raccolta dati
# FIX: usa in_basis=True così gli indici sono sempre 0,1,...,h11-1
# ─────────────────────────────────────────────

polytopes = fetch_polytopes(h11=2, as_list=True, limit=100)
print("Politopi trovati:", len(polytopes))

cy_data = {}

for i, P in enumerate(polytopes):
    try:
        t = P.triangulate(include_points_interior_to_facets=True)
        CY = t.get_cy()

        h11 = CY.h11()
        h21 = CY.h21()

        
        inter = CY.intersection_numbers(in_basis=True)

        # second_chern_class() restituisce i valori già nella base di h11
        # i primi h11 valori sono c2·J_i per i=0,...,h11-1
        c2 = CY.second_chern_class(in_basis=True)

        key = (h11, h21)
        if key not in cy_data:
            cy_data[key] = []

        cy_data[key].append((i, inter, c2))

    except Exception as e:
        print(f"  Politopo {i} skippato: {e}")
        continue

print("\nChiavi con più di una CY:")
for key in cy_data:
    if len(cy_data[key]) > 1:
        print(f"  {key}  →  {len(cy_data[key])} varietà")


####OUTPUT#####
Politopi trovati: 36

Chiavi con più di una CY:
  (2, 74)  →  2 varietà
  (2, 86)  →  9 varietà
  (2, 106)  →  3 varietà
  (2, 122)  →  3 varietà
  (2, 128)  →  5 varietà
  (2, 132)  →  2 varietà



------------------------------ CELLA 2-----------------------------------
def to_tensor(d, h11):
    """
    Converte il dizionario degli intersection numbers in un tensore NumPy 3D simmetrico.

    Con in_basis=True, le chiavi hanno sempre indici in {0,...,h11-1},
    quindi il tensore è esattamente di shape (h11, h11, h11).
    """
    t = np.zeros((h11, h11, h11), dtype=int)

    for (a, b, c), val in d.items():
        # Simmetrizza: assegna val a tutte le permutazioni di (a,b,c)
        for p in set(product((a, b, c), repeat=3)):
            if sorted(p) == sorted((a, b, c)):
                t[p] = val

    return t


def trova_matrice_equivalenza(cy1_data, cy2_data, h11=2, limit=5):
    """
    Cerca una matrice M ∈ GL(h11, ℤ) tale che:
      - M · c2_1 = c2_2
      - Σ_{xyz} M_ix M_jy M_kz κ1_xyz = κ2_ijk  per ogni i,j,k

    Se trovata restituisce M, altrimenti None.

    Parametri:
      limit: entrate di M cercate in [-limit, +limit]
             aumenta se non trova nulla (default 5 per sicurezza)
    """
    idx1, inter1, c2_1 = cy1_data
    idx2, inter2, c2_2 = cy2_data

    T1 = to_tensor(inter1, h11)
    T2 = to_tensor(inter2, h11)

    # Prendi i primi h11 valori di c2 (quelli nella base scelta)
    V1 = np.array(c2_1[:h11], dtype=int)
    V2 = np.array(c2_2[:h11], dtype=int)

    indices = range(h11)

    for entries in product(range(-limit, limit + 1), repeat=h11 * h11):
        M = np.array(entries, dtype=int).reshape(h11, h11)

        # Solo M ∈ GL(h11, ℤ): det = ±1
        if int(round(np.linalg.det(M))) not in [1, -1]:
            continue

        # Filtro rapido: verifica c2 prima 
        if not np.array_equal(M.dot(V1), V2):
            continue

        # Verifica trasformazione tensore di intersezione
        # κ'_ijk = Σ_{xyz} M_ix M_jy M_kz κ_xyz
        match = True
        for i, j, k in product(indices, repeat=3):
            val = sum(
                M[i, x] * M[j, y] * M[k, z] * T1[x, y, z]
                for x, y, z in product(indices, repeat=3)
            )
            if val != T2[i, j, k]:
                match = False
                break

        if match:
            return M

    return None


----------------------CELLA 3-------------------------------------
target_key = (2, 86)
lista_cands = cy_data[target_key]
h11 = target_key[0]

print(f"Analizzando {len(lista_cands)} varietà con Hodge {target_key}...\n")

found = False

for cy_a, cy_b in combinations(lista_cands, 2):
    M = trova_matrice_equivalenza(cy_a, cy_b, h11=h11, limit=10)

    if M is not None:
        print("-" * 50)
        print("TROVATE DUE VARIETÀ EQUIVALENTI!")
        print(f"Varietà A (Politopo {cy_a[0]})  <--->  Varietà B (Politopo {cy_b[0]})")
        print("\nMatrice di cambio base M:")
        print(M)
        print(f"\ndet(M) = {int(round(np.linalg.det(M)))}")
        print("\nVerifica c2:")
        V1 = np.array(cy_a[2][:h11])
        V2 = np.array(cy_b[2][:h11])
        print(f"  c2(A) = {V1}")
        print(f"  M·c2(A) = {M.dot(V1)}")
        print(f"  c2(B) = {V2}")
        print(f"  Coincidono: {np.array_equal(M.dot(V1), V2)}")
        print("-" * 50)
        found = True
        # Rimuovi il break se vuoi trovare TUTTE le coppie equivalenti

if not found:
    print("Nessuna coppia equivalente trovata.")
    print("Prova ad aumentare il parametro 'limit' in trova_matrice_equivalenza.")

##################OUTPUT##################################
Analizzando 9 varietà con Hodge (2, 86)...

--------------------------------------------------
TROVATE DUE VARIETÀ EQUIVALENTI!
Varietà A (Politopo 6)  <--->  Varietà B (Politopo 7)

Matrice di cambio base M:
[[2 1]
 [1 1]]

det(M) = 1

Verifica c2:
  c2(A) = [24  8]
  M·c2(A) = [56 32]
  c2(B) = [56 32]
  Coincidono: True