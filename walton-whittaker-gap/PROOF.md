# The window [9/8, 23/20] is a spectral gap of every cluster of the charged Walton–Whittaker tiling

## 0. Result

Let λ₁ = 9/8 and λ₂ = 23/20. For **both chains and every k ≥ 1** the matrix A(C_k) has exactly m = |C_{k-1}|
eigenvalues ≥ λ₁; all of them are > λ₂; none lies in [λ₁, λ₂].

* Spectrum of A(C_k): |C_{k-1}| eigenvalues above the window, |C_k| − |C_{k-1}| below.
* Gap label (value of the integrated density of states in the gap), finite k: ℓ_k = 1 − |C_{k-1}|/|C_k|
  * chain A (|C_k| = 1, 4, 16, 63, 245, 941, 3583, …): 3/4, 3/4, 47/63, 26/35, 696/941, 2642/3583, …
  * chain B (|C_k| = 2, 10, 42, 165, 632, 2396, 9049, …): 4/5, 16/21, 41/55, 467/632, 441/599, 6653/9049, …
* Both chains obey s_{k+6} = 8s_{k+5} − 22s_{k+4} + 25s_{k+3} − 7s_{k+2} − 5s_{k+1} + 2s_k, characteristic polynomial
  (x−1)(x⁵ − 7x⁴ + 15x³ − 10x² − 3x + 2), whose Perron root is θ = 3.772538098857… (the quintic is irreducible).
  Hence **ℓ = 1 − 1/θ = 0.734926467594…** below the gap and **1/θ = 0.265073532406…** above it, the roots of
  2ℓ⁵ − 7ℓ⁴ − 2ℓ³ + 13ℓ² − 9ℓ + 2 = 0 and 2y⁵ − 3y⁴ − 10y³ + 15y² − 7y + 1 = 0 respectively.

The proof follows the route of the hints: a lower-edge certificate from a matrix P (here P = 2 Σ u_v u_vᵀ), an upper-edge
certificate from p(x) = (2x²+4x+1)²(x − 23/20), exact rational arithmetic on the finite set of inflation cells,
and the min–max principle. The one point that needs care is the passage from cells to whole clusters (§4.4).
Everything is machine-checked in exact arithmetic by `python -m ww.verify` (≈ 35 s).

## 1. Setting

Sites are lattice points with a label (r, τ). The inflation σ sends a site p with label (r, τ) to four sites
2p, 2p+D_{r+3} (label (r,τ)), 2p+D_{r+2} (label (r−1,τ)) — the *triangle* — and 2p+D_{r+5} (label (r+2,−τ)).
σ commutes with the group G = Z₆ × Z₂ (rotation by 60° with r ↦ r+1; species flip). The tiling is Λ = ⋃_K σ^K({0}); other
seeds give G-images, so everything below is taken modulo G. The graph joins lattice neighbours of equal species;
C_k is a finite connected component. As in the statement, C_k consists of the m = |C_{k-1}| triangles T_v (v ∈ C_{k-1}) and of
*extras*: flipped fourth children w₃ of neighbouring sites w of the other species. "Off the boundary" is taken to mean
that every lattice neighbour of every site of the seed cluster C₀ is a site of the tiling; C₀ is one site (chain A) or an
adjacent pair (chain B). (The seeds that touch a hole of the tiling start further families, see §6.)
The only facts about the chains used below are the two from the statement: the triangles of C_k are exactly the T_v, v ∈ C_{k-1}, and every other
vertex of C_k is a flipped fourth child of a lattice neighbour of some site of C_{k-1}.

Notation: n = |C_k|, A = A(C_k), u_v ∈ ℝⁿ the indicator of T_v, U = [u_v], S = span U (dim S = m).

## 2. Two counting lemmas

**Lemma 1.** Let A be real symmetric, λ > 0, P real symmetric with Q := λI − A + APA ≻ 0. Then A has at most n₊(P)
eigenvalues ≥ λ.

*Proof.* Let E be the span of the eigenvectors with eigenvalue ≥ λ, e = dim E, and N the span of the eigenvectors of P
with eigenvalue ≤ 0, dim N = n − n₊(P). A is injective on E (λ > 0), so dim AE = e. If e > n₊(P) then
dim AE + dim N > n and there is x ∈ E, x ≠ 0, with y = Ax ∈ N. Then
xᵀQx = (λ|x|² − xᵀAx) + yᵀPy ≤ 0 + 0, contradicting Q ≻ 0. ∎

**Lemma 2.** If p(t) = q(t)²(t − λ₂) with q real and S is a subspace with xᵀp(A)x > 0 for all x ∈ S∖{0}, then A has at
least dim S eigenvalues > λ₂.

*Proof.* Let F be the span of the eigenvectors with eigenvalue ≤ λ₂. On F, xᵀp(A)x = Σ p(μ_i) ξ_i² ≤ 0 because
p(μ) ≤ 0 for μ ≤ λ₂. So S ∩ F = 0, dim F ≤ n − dim S. ∎

Take P = 2 Σ_v u_v u_vᵀ (positive semidefinite, n₊(P) = m; P − β(I − Π_S) with β = 10⁻³ is invertible and works equally,
because ‖A‖ ≤ 6) and q(x) = 2x² + 4x + 1. Put

* **Q₁ = (9/8)I − A + APA = (9/8)I − A + 2 Σ_v (A u_v)(A u_v)ᵀ,**
* **M = Uᵀ p(A) U,  p(x) = (2x²+4x+1)²(x − 23/20)** (so that cᵀMc = xᵀp(A)x for x = Σ c_v u_v).

**Proposition.** For both chains and all k ≥ 1: Q₁ ⪰ 0.0799·I and M ⪰ 4.148·I.

*The Proposition implies the Result:* Lemma 1 (λ = 9/8) gives #{μ ≥ λ₁} ≤ m; Lemma 2 (S = span U) gives
#{μ > λ₂} ≥ m. Hence m ≤ #{μ > λ₂} ≤ #{μ ≥ λ₁} ≤ m: all inequalities are equalities, so there is no eigenvalue in [λ₁, λ₂]. ∎

## 3. Why these objects

On S the triangle has eigenvalue 2, so S carries the upper band, but the Rayleigh quotients of A on S reach down to
≈ 0.29: the true band eigenvectors are not in S. Lemma 2 repairs this by testing p(A) instead of A; the double factor q²
suppresses the bulk components of the triangle vectors (none of the degree-1 candidates g = x + c, c ∈ {0, …, 3}, works already on
the cluster C₃ of chain B with 165 vertices; the degree-2 factor q of the hint does). Lemma 1 with P ≈ Π_S is the statement
"A < 9/8 on {x : u_vᵀAx = 0 ∀v}".

## 4. Proof of the Proposition: finite verification

### 4.1 Cells

For a site z let W(z) be its *window*: the labels of all sites within lattice distance 7 of z (positions relative to z).
Two facts, checked by `ww.core.geometric_check` and `tests/equivariance.py`: (i) the radius-7 window of a child of z is a
function of W(z) (a site within distance 7 of a child has its parent within distance 7 of z), (ii) this map commutes with G.
One inflation step applied to W(z) therefore determines all sites of the next level within Euclidean distance 10 of 2z, in particular
the triangle T_z, the nearest-neighbour graph around it, the roles (child index, parent) of all these sites, and the windows
of the possible new members.

*All windows.* Iterating the children map from the window of the origin gives every window of Λ: 3542 windows mod G.
*Seeds.* The windows centred at a singleton cluster / at a site of a two-site cluster all of whose lattice neighbours are sites:
103 (chain A) + 22 (chain B).
*Cells.* The least G-invariant set 𝕎 of windows containing the seeds and closed under the map

  W(z) ↦ { windows of the three vertices of T_z, and of the flipped fourth children w₃ of **all** neighbours w of z of the other species }

has **3223** elements (stable after 8 rounds); with the seeds adjacent to holes of the tiling included it has 3440.
(Taking *all* fourth children, not only those that turn out to lie in the cluster of T_z, keeps the cell set independent of any
decision about which sites belong to a given cluster, which for large clusters is not determined inside a window.)

**Lemma 3 (coverage).** For both chains, every k ≥ 0 and every z ∈ C_k: W(z) ∈ 𝕎 (mod G).

*Proof.* Induction on k. k = 0: seeds. Step: a vertex of C_{k+1} is a vertex of T_z for some z ∈ C_k, or a flipped fourth child w₃ of a
neighbour w of some z ∈ C_k (structure of the chains); its window is produced from W(z) by the map above. ∎

### 4.2 Local sum-of-squares certificates

*Lower edge.* For a vertex c of C_k let B(c) be the vertices at graph distance ≤ 2 from c; its type is
τ(c) = (child index of c, G-orbit of the radius-2 window of the parent of c). There are 896 types; a type fixes the shape of B(c) in a
reference frame (checked). To each type belongs a positive definite rational matrix G_τ on B_τ (sizes 4…15). Put
Q̂ = Σ_c ι_cᵀ G_{τ(c)} ι_c and R = Q₁ − Q̂.

*Upper edge.* For a coarse site v ∈ C_{k-1} (↔ triangle T_v) let B_H(v) be the triangles within distance ≤ 2 of T_v in the
triangle-adjacency graph H of C_k; type τ(v) = G-orbit of the radius-4 window of v. 1028 types, positive definite rational
G_τ (sizes 1…13), M̂ = Σ_v ι_vᵀ G_{τ(v)} ι_v, R = 20M − M̂.

### 4.3 What is verified in exact arithmetic (`python -m ww.verify`)

1. every G_τ is positive definite (exact elimination over ℚ, all pivots > 0);
2. for every cell in 𝕎 and every row x it owns (Q₁: the 3 triangle vertices of z and the fourth children of all other-species
   neighbours of z, 18473 rows in total; M: the row of z, 3223 rows) the row of R satisfies the Gershgorin inequality
   R_xx − Σ_{y≠x} |R_xy| ≥ ρ, with

   * ρ_L = 85844911 / 1073741824 = **0.07995** for Q₁,
   * ρ_U = 696063391 / 8388608 = **82.977** for 20·M, i.e. 4.1489 for M.

No floating point enters this step. (`--all-seeds` repeats it on the 3440 cells generated by all seeds, boundary ones included.)
The row of a vertex x of a cluster is the row of the same-species graph around x (the ball of radius 3 does not depend on which
cluster x is known to belong to), so rows of fourth children are checked whether or not they lie in the cluster of z.

### 4.4 Lifting to whole clusters

Let C_k be any cluster of either chain, F ∈ {Q₁, 20M} its form, and F = F̂ + R with F̂ = Σ ι_cᵀ G_{τ(c)} ι_c.
Each summand of F̂ is positive semidefinite (item 1), so F̂ ⪰ 0. R is symmetric and, by Lemma 3 and item 2, every row of R of the cluster is
a row of some cell, hence diagonally dominant with margin ρ, so R ⪰ ρI by Gershgorin. Therefore F ⪰ ρ I, i.e.
Q₁ ⪰ 0.0799 I and M ⪰ 4.1489 I. ∎

(Remark: positive definiteness of the cell blocks of F alone would *not* lift — interlacing gives the opposite inequality.
What lifts is a decomposition into positive semidefinite local pieces with a locally dominant remainder; the pieces were found
by a semidefinite program, which is irrelevant for the proof — only the exact check of §4.3 matters.)

## 5. Gap labels

Let N_k(τ) be the number of sites of C_k whose cell is τ. A site of C_{k-1} contributes its 3 triangle vertices; each extra is
attributed to one adjacent site of C_{k-1} of the same cluster (the one minimising (dir(w→z′) − r_w) mod 6). The cell of z
determines the cells of everything attributed to it, so N_{k+1} = 𝕄 N_k with an explicit non-negative integer matrix 𝕄 on the cells
(`ww/counts.py`). From a chain-A seed the Krylov space of 𝕄 has exact dimension 12, so |C_k| = 1ᵀ𝕄ᵏN₀ satisfies a linear
recurrence for all k; the minimal one has order 6 (checked to k = 60 for both chains; all seeds of a chain give the same sequence).
Its dominant root θ gives lim |C_{k-1}|/|C_k| = 1/θ, the fraction of the spectrum above the window; ℓ = 1 − 1/θ is the label of the gap.
(Here "same cluster" is decided inside the window; the resulting sizes agree with a brute-force construction of the clusters for
k ≤ 9. This bookkeeping is used only for the label, not for the gap.)

## 6. Reproduction, independent checks, scope

```
python -m ww.verify              # the proof: closure of cells, exact check of both certificates (≈ 35 s)
python -m ww.verify --all-seeds  # same on all 3440 cells (also chains starting at holes of the tiling)
python -m ww.counts              # exact |C_k|, recurrence, gap labels
python -m tests.crosscheck 6 5   # brute-force clusters: eigenvalue counts, global exact residuals (float = not the proof)
python -m tools.gen_constraints L 2 L.pkl --all-seeds; python -m tools.solve_sdp L.pkl cert_L.json.gz L   # regenerate certificates (cvxpy, scs)
```

Numerics on actual clusters (dense eigensolver, k ≤ 6 for A and k ≤ 5 for B): #eig ≥ 9/8 = |C_{k-1}|, none in [9/8, 23/20];
largest eigenvalue below the window ≤ 1.0071, smallest above it 1.2549 (B, k = 5, slowly decreasing). The window is therefore
a sub-window of a larger gap; the certificates are not optimal and the margins are small, as the hints said.

*Scope.* Trusted inputs: the two structural facts quoted in §1 and the correctness of the code in `ww/` (cross-checked against an
independent brute-force construction of the actual clusters up to level 8: all member windows lie in 𝕎, cell rows equal the true rows,
global exact residual margins reproduce the certified ones). The seeds adjacent to holes of the tiling (61 further singleton
windows) start three more families of finite clusters, with sizes (1, 4, 16, 62, 236, …), (1, 4, 15, 56, 210, …), (1, 4, 15, 55, 201, …);
they are "on the boundary" and not part of the statement, but the same certificates verify on their cells too.
