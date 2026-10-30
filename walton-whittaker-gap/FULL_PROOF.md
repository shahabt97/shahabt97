# A spectral gap [9/8, 23/20] for the clusters of the charged Walton–Whittaker tiling

**A complete, self-contained, computer-assisted proof (exact rational arithmetic)**

This file contains everything needed to evaluate the proof: the problem and its precise reading (§1), the statements (§2), every
mathematical argument (§3–§7), an honest list of assumptions and risks (§9), and — verbatim — all source code, all program
outputs and fingerprints of the data files (appendices). Nothing outside this file is needed to *read* the proof; the finite
computations can be re-run with the code in Appendix A (Python ≥ 3.9, numpy).

---

## Contents

0. Summary, architecture, claim ledger
1. The problem and its formalization
2. Results (Theorem A: the gap; Theorem B: gap labels)
3. Paper layer I — counting and lifting lemmas; reduction to two positivity statements
4. Paper layer II — symmetry, windows, finiteness, locality, and the certificate lemma
5. The exact computations
6. Validation (not used in the proof)
7. Gap labels (proof of Theorem B)
8. Generality (chains that start next to holes)
9. Assessment: assumptions, risks, what is *not* proved
- Appendix A — source code (verbatim)
- Appendix B — program outputs (verbatim)
- Appendix C — fingerprints of code and data
- Appendix D — how to reproduce

---

## 0. Summary, architecture, claim ledger

### 0.1 What is proved

Let **λ₁ = 9/8** and **λ₂ = 23/20**. For **both chains** of finite clusters of the charged Walton–Whittaker tiling and **every k ≥ 1**
the adjacency matrix A(C_k) has exactly m = |C_{k−1}| eigenvalues ≥ λ₁; all of them are > λ₂; none lies in [λ₁, λ₂]
(Theorem A). The integrated-density-of-states value in the gap, i.e. the gap label, is 1 − |C_{k−1}|/|C_k| at level k and
**1 − 1/θ = 0.7349264675940…** in the limit (above the gap: 1/θ = 0.2650735324059…), where θ = 3.7725380988570… is the Perron root of
x⁵ − 7x⁴ + 15x³ − 10x² − 3x + 2 (Theorem B).

### 0.2 Architecture

1. **Two counting lemmas** (§3.1–3.2, pure linear algebra, short proofs).
   *Lemma 1* (inertia): if (9/8)I − A + APA ≻ 0 for some symmetric P, then A has at most n₊(P) eigenvalues ≥ 9/8.
   *Lemma 2* (min–max): if xᵀp(A)x > 0 on an m-dimensional subspace S, where p(t) = q(t)²(t − 23/20), then A has at least m eigenvalues > 23/20.
   With P = 2 Σ_v u_v u_vᵀ (u_v = indicator of the triangle of the site v ∈ C_{k−1}; n₊(P) = m) and S = span{u_v}, the two lemmas squeeze
   the number of eigenvalues above the window to exactly m and exclude the window.
2. **Two positivity statements** (Proposition P, §3.4) for **all k and both chains at once**:
   Q₁ := (9/8)I − A + 2Σ_v (Au_v)(Au_v)ᵀ ⪰ 0.0799·I and M := Uᵀp(A)U ⪰ 4.148·I, with p(x) = (2x²+4x+1)²(x − 23/20).
3. **Local-to-global lifting** (Lemma 3): a symmetric matrix F that equals a sum of *positive semidefinite local pieces*
   plus a *diagonally dominant* remainder is positive definite. (Positivity of cell blocks alone would not lift.)
4. **Finiteness** (§4): the local structure of every cluster of both chains, to the radius needed, is described by a finite set 𝕎 of
   "inflation cells" (radius-7 neighbourhoods modulo the symmetry group of the tiling), closed under inflation; 𝕎 has 3223 elements.
5. **Exact certificate** (§4.6, §5): 896 + 1028 rational positive definite matrices ("pieces") and, for each of the 3223 cells, the row inequalities
   of the remainder, checked in exact rational arithmetic. Result: ρ_L = 85844911/1073741824 ≈ 0.07995 and ρ_U/20 = 696063391/167772160 ≈ 4.1489.

### 0.3 Status tags

* **[H]** hypothesis taken from the problem statement (the structure of the chains, §1.3);
* **[P]** proved on paper in this document;
* **[X]** exact computation (python integers / `fractions.Fraction`; code in Appendix A, output in Appendix B);
* **[V]** validation by brute force / floating point; **not used** in the proof.

| # | Claim | Tag | Where |
|---|-------|-----|-------|
| H0–H2 | structure of the chains (clusters; triangles of C_k = triangles of the sites of C_{k−1}; every other vertex is a flipped fourth child of a neighbour of a site of C_{k−1}) | [H] | §1.3 |
| L1, L2, L3 | counting lemmas, lifting lemma | [P] | §3 |
| L4 | inflation commutes with the symmetry group G | [P] (+[V] test) | §4.1 |
| L5 | the window of a child is a function of the window of its parent; determined region | [P] | §4.2 |
| L6 | the tiling is collision free; unique parents | [P] using the closure computation [X] | §4.3 |
| L7 | all windows of the tiling lie in a finite computed set 𝒲 (3542 classes) | [P] + [X] | §4.4 |
| L8 | every site of every chain cluster has its window in the finite set 𝕎 of cells (3223 classes) | [P] given [H] + [X] | §4.5 |
| L9 | locality: every quantity in a row inequality is a function of one cell | [P] (distance bounds) + [X] guards | §4.6 |
| L10 | certificate ⇒ Proposition P | [P] | §4.7 |
| C1 | computation of 𝒲, seeds, 𝕎 | [X] | §5 |
| C2 | all 896 + 1028 pieces are positive definite | [X] | §5 |
| C3 | all row inequalities hold on all cells | [X] | §5 |
| T-A | Theorem A | follows from L1, L2, P | §3.5 |
| T-B | gap labels | finite k: [P] from T-A; limit: [P] given the recurrence; recurrence: [X] in the substitution model (+[V]); algebra of the quintic: [X] exact (`tools/algebra.py`) | §7 |

### 0.4 Where the weight of the argument lies (for a judge)

* Lemmas 1–3 (§3), the reduction to Proposition P (§3.5), Lemmas 4–6 and 8–10 (§4) and the algebra of §7 are short proofs that can be checked by hand.
* The computer-assisted part is a *finite* statement (C1–C3, §5) whose meaning is spelled out in §4.5–§4.7 and §5.3 and whose code is in Appendix A. A suggested order for checking it: (1) the definitions of the two local forms in `ww/cells.py` (A.3) against §4.6; (2) the exact checker `ww/verify.py` (A.4) against Lemma 10; (3) the closure `ww/closure.py` (A.2) against §4.4–§4.5; (4) the tiling rules `ww/core.py` (A.1) against §1.2.
* The main logical risks are in §9: the reading of the problem (H0–H2), the correctness of the code, and — only for the limit in Theorem B — the substitution model of §7.
* Independent brute-force validation of the code (a different code path, in floating point and in exact arithmetic) is reported in §6 and Appendix B; it is **not** part of the proof.

---

## 1. The problem and its formalization

### 1.1 The problem as posed

> The charged Walton–Whittaker tiling is on the triangular lattice. Each site is a label (r,τ) with r∈ℤ/6ℤ and species τ=±1.
> Under charged midpoint inflation, a site v labelled (r,τ) has four children:
> • 2v and 2v+D_{r+3}, both labelled (r,τ)
> • 2v+D_{r+2}, labelled (r+5,τ)
> • 2v+D_{r+5}, labelled (r+2,−τ)
> D_0,…,D_5 are the six lattice directions in counterclockwise order. The first three children form a triangle. Join nearest
> neighbours of the same species. Every finite cluster off the boundary lies on one of two chains. A level-k cluster C_k, on either
> chain, is the triangles of all sites of C_{k-1}, plus the flipped fourth children of neighbouring sites of the other species,
> attached as a matching. A(C_k) is the adjacency matrix of that cluster.
> Prove: there exist λ_1<λ_2 such that for every k≥1 and both chains, A(C_k) has exactly |C_{k-1}| eigenvalues above λ_1 and none in
> [λ_1,λ_2]. A float is not a proof. Give the exact interval and the gap labels it implies.

Hints that came with the problem: the window is [9/8, 23/20]; for the lower edge find an invertible P such that (9/8)I − A + APA is positive
definite on each inflation cell; for the upper edge take a polynomial with a simple factor (x − 23/20), of degree 5, a constant times
(2x²+4x+1)²(x − 23/20), and show that its quadratic form is positive definite on the complement of the "bulk subspace"; check both forms in exact
rational arithmetic on the finite set of inflation cells and lift with the min–max principle.

### 1.2 Conventions

**Lattice.** A lattice point is a pair (a,b) ∈ ℤ², standing for a·D₀ + b·D₁ with D_j = e^{iπj/3}. Thus
D₀ = (1,0), D₁ = (0,1), D₂ = (−1,1), D₃ = (−1,0), D₄ = (0,−1), D₅ = (1,−1); indices are mod 6, D_{j+3} = −D_j and D_{j+1} = D_j + D_{j+2}.
hexd(a,b) := max(|a|, |b|, |a+b|) is the lattice (graph) distance in the triangular lattice; it is a norm on ℤ² (so hexd(2q) = 2·hexd(q) and
the triangle inequality holds), and its unit sphere is {±D_j}. The six *lattice neighbours* of p are p + D_j. Rotation by 60° is ρ(a,b) = (−b, a+b) (ρD_j = D_{j+1}).

**Sites, inflation.** A *site* is a lattice point p with a label ℓ(p) = (r,τ) ∈ ℤ/6 × {±1}. For a site p with label (r,τ) put
ε₀ = 0, ε₁ = D_{r+3}, ε₂ = D_{r+2}, ε₃ = D_{r+5} and λ₀ = λ₁ = (r,τ), λ₂ = (r+5,τ), λ₃ = (r+2,−τ). The *children* of p are
c_j(p) := 2p + ε_j with label λ_j (j = 0,1,2,3). The *triangle* of p is T_p := {c₀(p), c₁(p), c₂(p)}; its three points are pairwise lattice neighbours
(the differences are D_{r+3}, D_{r+2}, D_{r+2} − D_{r+3} = D_{r+1}) and have the same species τ. The child c₃(p) has species −τ ("flipped").
The inflation of a labelled set S is σ(S) := {children of the sites of S}.

**The tiling.** Λ_K := σ^K({0}), where 0 carries the label (0,+1); Λ := ⋃_K Λ_K. (Other seed labels give images of Λ under the group G below;
nothing changes.) Lemma 6 shows that Λ is a well-defined labelled set with σ(Λ) = Λ and that every y ∈ Λ has a unique *parent* π(y) ∈ Λ and *child index*
j(y) ∈ {0,1,2,3} with y = c_{j(y)}(π(y)) (the origin is its own parent).

**Symmetry group.** G = ℤ₆ × ℤ₂ acts on labelled lattice points: the generator of ℤ₆ maps (p,(r,τ)) ↦ (ρp,(r+1,τ)); the generator of ℤ₂ maps (p,(r,τ)) ↦ (p,(r,−τ)).

**The graph and clusters.** 𝒢 is the graph on Λ whose edges join lattice neighbours of equal species ("nearest neighbours of the same species").
A *cluster* is a connected component of 𝒢. For a set C ⊂ Λ, A(C) is the adjacency matrix of the induced subgraph 𝒢[C].
Graph distance in 𝒢 is ≥ hexd.

**Windows.** For a site p and R ≥ 1, the *window* W_R(p) is the labelled set {(q, ℓ(p+q)) : hexd(q) ≤ R, p+q ∈ Λ} (positions relative to p; absent lattice points are recorded as absent).
G acts on windows (rotating the offsets and recoding the labels); windows are always considered *modulo G* ("G-classes").

### 1.3 The two chains: hypotheses H0–H2

The statement says that "every finite cluster off the boundary lies on one of two chains" and that C_k is "the triangles of all sites of C_{k−1},
plus the flipped fourth children of neighbouring sites of the other species". I read this as follows; **these are the only facts about the chains that
the proof uses**, and they are taken from the statement [H].

* **(H0)** C₀ is either a single site {z} with no same-species lattice neighbour (**chain A**), or two lattice-adjacent sites {z, z′} of equal species, neither having another
  same-species neighbour (**chain B**); in both cases *every lattice neighbour of every site of C₀ is a site of Λ* ("off the boundary").
  For k ≥ 1, C_k is a finite cluster.
* **(H1)** For k ≥ 1: {y ∈ C_k : j(y) ∈ {0,1,2}} = ⋃_{v ∈ C_{k−1}} T_v. (The triangle vertices of C_k are exactly the triangles of the sites of C_{k−1}.)
* **(H2)** For k ≥ 1: every y ∈ C_k with j(y) = 3 (an *extra*) is c₃(w) for some lattice neighbour w of some v ∈ C_{k−1}. (Then w has the species opposite to that of C_k, i.e. "the other species".)

Consequences used repeatedly: the T_v (v ∈ C_{k−1}) are disjoint triangles, so C_k contains m = |C_{k−1}| triangles and n = 3m + (number of extras) vertices.

*Why this is the right reading.* A census of the patch Λ₉ (output in Appendix B, §B.6) finds exactly the sizes
{1, 2, 4, 10, 16, 42, 63, 165, 245, 941, 3583} among the finite clusters all of whose lattice neighbours are sites; these are the first terms of the two size
sequences (1, 4, 16, 63, …) and (2, 10, 42, 165, …) of chains A and B, and nothing else occurs [V]. H1 and H2 were also checked by brute force on the actual clusters
up to k = 9 (chain A) and k = 8 (chain B) (Appendix B, §B.7) [V]. The phrase "attached as a matching" is not used: in the actual clusters an extra is adjacent to 0, 1, 2 or 3 triangle vertices and to at most one other extra (Appendix B, §B.12).

### 1.4 Remarks on conventions

* *Conventions.* The only choices are the identification D_j = e^{iπj/3} with D₀ = (1,0) and the label (0,+1) of the origin. Starting the labelling of directions at another lattice direction amounts to adding a constant to all r, i.e. to a relabelling by an element of G (a rotation), which changes nothing below.
* "Nearest neighbours" = lattice neighbours (distance 1). "Above λ₁" and "≥ λ₁" are equivalent in Theorem A because no eigenvalue lies in [λ₁, λ₂].
* "|C_{k−1}|" is the number of sites of C_{k−1}.

---

## 2. Results

**Theorem A (the gap).** Let λ₁ = 9/8, λ₂ = 23/20. Under H0–H2, for both chains and every k ≥ 1, the matrix A(C_k) has exactly m = |C_{k−1}| eigenvalues ≥ λ₁,
all of them are > λ₂, and no eigenvalue lies in [λ₁, λ₂]. Equivalently, spec A(C_k) ⊂ (−∞, λ₁) ∪ (λ₂, ∞) with |C_k| − |C_{k−1}| eigenvalues below and |C_{k−1}| above.

**Theorem B (gap labels).** For E ∈ [λ₁, λ₂] the normalised counting function N_k(E) = #{eigenvalues of A(C_k) below E}/|C_k| equals 1 − |C_{k−1}|/|C_k|.
Both chains satisfy s_{k+6} = 8s_{k+5} − 22s_{k+4} + 25s_{k+3} − 7s_{k+2} − 5s_{k+1} + 2s_k, whose characteristic polynomial is (x−1)(x⁵ − 7x⁴ + 15x³ − 10x² − 3x + 2);
hence N_k(E) → ℓ := 1 − 1/θ = 0.734926467594067…, with θ = 3.772538098857050… the Perron root of the (irreducible) quintic. Thus ℓ is the root of
2ℓ⁵ − 7ℓ⁴ − 2ℓ³ + 13ℓ² − 9ℓ + 2 = 0 near 0.7349, and the fraction of the spectrum above the gap is 1/θ = 0.265073532405…, the root of 2y⁵ − 3y⁴ − 10y³ + 15y² − 7y + 1 = 0 near 0.265.
Finite-k labels: chain A: 3/4, 3/4, 47/63, 26/35, 696/941, 2642/3583, 9991/13574, …; chain B: 4/5, 16/21, 41/55, 467/632, 441/599, 6653/9049, 25086/34135, ….

(Theorem A is proved in §3–§5 under H0–H2 and the computer-assisted steps C1–C3. In Theorem B the finite-k labels follow directly from Theorem A, while the limit also uses the size recurrence; §7 says exactly what is proved and what is only validated.)

---

## 3. Paper layer I — counting and lifting lemmas

Throughout, A is a real symmetric n×n matrix, n₊(P) is the number of positive eigenvalues of a symmetric matrix P, and ⪰/≻ are the Loewner orders.

### 3.1 Lemma 1 (inertia bound)

**Lemma 1.** Let λ > 0 and let P be real symmetric with Q := λI − A + APA ≻ 0. Then A has at most n₊(P) eigenvalues ≥ λ.

*Proof.* Let E be the span of the eigenvectors of A with eigenvalue ≥ λ, e = dim E, and let N be the span of the eigenvectors of P with eigenvalue ≤ 0, so dim N = n − n₊(P).
Since λ > 0, A is injective on E, hence dim AE = e. Suppose e > n₊(P). Then dim AE + dim N > n, so there is x ∈ E, x ≠ 0, with y := Ax ∈ N. Then
xᵀQx = (λ|x|² − xᵀAx) + yᵀPy. The first bracket is ≤ 0 because x ∈ E, the second is ≤ 0 because y ∈ N. So xᵀQx ≤ 0, contradicting Q ≻ 0. ∎

(Invertibility of P is not needed. If one insists on an invertible P, replace P by P_β = P − β(I − Π) with Π the orthogonal projector onto range P, β = 10⁻³: n₊ is unchanged,
and Q changes by −β A(I − Π)A ⪰ −β‖A‖² I with ‖A‖ ≤ 6 (maximal degree of the triangular lattice), so the margin 0.0799 of Proposition P below survives: 36β = 0.036 < 0.0799.)

### 3.2 Lemma 2 (min–max bound)

**Lemma 2.** Let λ₂ ∈ ℝ, let q be a real polynomial, p(t) := q(t)²(t − λ₂), and let S be a subspace with xᵀp(A)x > 0 for all x ∈ S∖{0}. Then A has at least dim S eigenvalues > λ₂.

*Proof.* Let F be the span of the eigenvectors of A with eigenvalue ≤ λ₂. For x = Σ ξ_i φ_i ∈ F (orthonormal eigenvectors, eigenvalues μ_i ≤ λ₂): xᵀp(A)x = Σ q(μ_i)²(μ_i − λ₂) ξ_i² ≤ 0.
Hence S ∩ F = {0}, so dim S + dim F ≤ n, i.e. dim F ≤ n − dim S, i.e. at least dim S eigenvalues are > λ₂. ∎

### 3.3 Lemma 3 (lifting: positive semidefinite pieces + diagonally dominant remainder)

**Lemma 3.** Let V be a finite set and F a real symmetric matrix indexed by V. Suppose that for each c ∈ V there are a subset B_c ⊂ V and a real symmetric *positive semidefinite*
matrix G_c indexed by B_c, and put F̂ := Σ_c ι_cᵀ G_c ι_c (ι_c: ℝ^V → ℝ^{B_c} the coordinate restriction), R := F − F̂. If for every x ∈ V
R_xx − Σ_{y≠x} |R_xy| ≥ ρ, then F ⪰ ρ I.

*Proof.* For z ∈ ℝ^V: zᵀFz = Σ_c (ι_c z)ᵀ G_c (ι_c z) + zᵀRz ≥ zᵀRz. Since R is symmetric and 2|z_x z_y| ≤ z_x² + z_y²,
zᵀRz ≥ Σ_x R_xx z_x² − Σ_{x≠y} |R_xy| |z_x||z_y| ≥ Σ_x (R_xx − Σ_{y≠x}|R_xy|) z_x² ≥ ρ|z|². ∎

(Remark. Positive definiteness of principal blocks F_{BB} would *not* suffice — the matrix [[1,2],[2,1]] has positive definite 1×1 blocks and is indefinite. What lifts is a decomposition of F into positive semidefinite local pieces plus a locally diagonally dominant remainder.)

### 3.4 The two forms and Proposition P

Fix a chain, k ≥ 1, put A = A(C_k) (n = |C_k|), m = |C_{k−1}|, and for v ∈ C_{k−1} let u_v ∈ ℝ^{C_k} be the indicator vector of the triangle T_v ⊂ C_k (H1). The u_v have disjoint supports, so U = [u_v] has rank m and S := span{u_v} has dimension m;
the hint's "bulk subspace" is S^⊥ (vectors with zero sum on every triangle, extras free), and the form M below lives on its complement S. Define

* P := 2 Σ_v u_v u_vᵀ (positive semidefinite; its nonzero eigenvalues are 2|u_v|² = 6, multiplicity m; so n₊(P) = m);
* **Q₁ := (9/8)I − A + APA = (9/8)I − A + 2 Σ_v (Au_v)(Au_v)ᵀ**;
* q(x) := 2x² + 4x + 1, **p(x) := q(x)²(x − 23/20)**; since q² = 4x⁴ + 16x³ + 20x² + 8x + 1 and 20(x − 23/20) = 20x − 23,
  **20·p(x) = 80x⁵ + 228x⁴ + 32x³ − 300x² − 164x − 23** (integer coefficients);
* **M := Uᵀ p(A) U** (an m×m matrix indexed by C_{k−1}); for x = Σ_v c_v u_v one has xᵀp(A)x = cᵀMc.

**Proposition P.** For both chains and every k ≥ 1:
 (L) Q₁ ⪰ ρ_L·I with ρ_L = 85844911/1073741824 = 0.0799…;
 (U) M ⪰ ρ_U′·I with ρ_U′ = 696063391/167772160 = 4.1488….

Proposition P is proved in §4 (paper) and §5 (exact computation).

### 3.5 Proposition P ⇒ Theorem A

By (L), Q₁ ≻ 0, so Lemma 1 (λ = 9/8, P as above) gives #{eigenvalues of A ≥ 9/8} ≤ n₊(P) = m.
By (U), xᵀp(A)x = cᵀMc > 0 for all x = Uc ≠ 0 (U has full column rank), so Lemma 2 (S = span U, dim S = m) gives #{eigenvalues > 23/20} ≥ m. Therefore
m ≤ #{μ > 23/20} ≤ #{μ ≥ 9/8} ≤ m. All three are equal to m; since #{μ > 23/20} = #{μ ≥ 9/8}, no eigenvalue lies in [9/8, 23/20]. ∎

### 3.6 Remarks on the two forms (not used in the proof)

* Q₁ ≻ 0 says in particular that xᵀAx < (9/8)|x|² for every x ≠ 0 with u_vᵀAx = 0 for all v (for such x, APAx = 0 and xᵀQ₁x = (9/8)|x|² − xᵀAx).
* The polynomial p of the hint is q(x)²(x − 23/20) with q(x) = 2x² + 4x + 1; as in the hint, the squared factor q² is there to remove the bulk and the linear factor changes sign at the upper edge 23/20. Any positive lower bound for M = Uᵀp(A)U suffices for Lemma 2.

---

## 4. Paper layer II — symmetry, windows, finiteness, locality, certificate

### 4.1 Symmetry

**Lemma 4.** σ commutes with G; consequently G maps Λ_K, Λ, parents, child indices, triangles, windows and the graph 𝒢 to the corresponding objects of the image, and A(·), Q₁, M are invariant.

*Proof.* Rotation: for a site (p,(r,τ)) with image (ρp,(r+1,τ)): c₁(ρp) = 2ρp + D_{(r+1)+3} = ρ(2p + D_{r+3}) = ρ c₁(p); likewise c₂, c₃ (D_{(r+1)+2}, D_{(r+1)+5}); the labels: λ₀,λ₁: (r+1,τ) = image of (r,τ);
λ₂: ((r+1)+5, τ) = (r+6,τ) = image of (r+5,τ); λ₃: ((r+1)+2, −τ) = (r+3,−τ) = image of (r+2,−τ). Flip: obvious. ρ preserves lattice neighbours, both generators preserve "equal species". ∎

(An independent computer test of the window version of this statement is in `tests/equivariance.py`, Appendix B §B.5 [V].)

### 4.2 Windows of children (Lemma 5)

For a labelled set S, a site p ∈ S and R ≥ 1 let B_R(p) := {q ∈ ℤ² : hexd(q − p) ≤ R} and let P_R(p; S) be the pattern S ∩ B_R(p) (labels kept, positions taken relative to p); the window W_R(p) of §1.2 is P_R(p; Λ). For a pattern P centred at p let ν_j(P) be the pattern in B_R(c_j(p)) formed by the children of the sites of P (j = 0,…,3).

**Lemma 5.** Let S be a labelled set, p ∈ S, y := c_j(p), and y′ = c_i(u′) ∈ σ(S) with u′ ∈ S.
 (a) If hexd(y′ − y) ≤ R then hexd(u′ − p) ≤ (R + 2)/2.
 (b) Hence, for R ≥ 2: P_R(y; σ(S)) = ν_j(P_R(p; S)). The window of a child is a function of the window of its parent (and of j).
 (c) If hexd(y′ − 2p) ≤ 13 then hexd(u′ − p) ≤ 7. Hence σ(S) ∩ B₁₃(2p) equals the set of children of the sites of P₇(p; S), intersected with B₁₃(2p) (the *determined region* of p). In particular absence of a site in this region is determined by P₇(p; S).
 (d) If w ∈ S is a lattice neighbour of p and y := c₃(w), then hexd(y − 2p) ≤ 3 and P₇(y; σ(S)) is determined by P₇(p; S).
 (e) |q| ≥ (√3/2)·hexd(q) for q ∈ ℤ² (|·| is the Euclidean length for lattice spacing 1). In particular |q| ≤ 10 implies hexd(q) ≤ 11 ≤ 13.

*Proof.* Write 2u′ − y′ = −ε and y − 2p = ε′ with ε, ε′ ∈ {0} ∪ {unit vectors}, so hexd(ε), hexd(ε′) ≤ 1; recall that hexd is a norm and hexd(2q) = 2·hexd(q).
(a) 2(u′ − p) = −ε + (y′ − y) + ε′, so 2·hexd(u′ − p) ≤ R + 2.
(b) For R ≥ 2, (R + 2)/2 ≤ R, so every site of σ(S) in B_R(y) is a child of a site of P_R(p; S); conversely children of sites of S that lie in B_R(y) are sites of σ(S).
(c) 2(u′ − p) = −ε + (y′ − 2p), so 2·hexd(u′ − p) ≤ 1 + 13 = 14.
(d) y − 2p = 2(w − p) + ε₃(w), so hexd(y − 2p) ≤ 2 + 1. A site of σ(S) within distance 7 of y lies within distance 10 ≤ 13 of 2p; apply (c).
(e) a² + ab + b² = ¾a² + (a/2 + b)² ≥ ¾a²; symmetrically ≥ ¾b², and a² + ab + b² = ¾(a+b)² + ¼(a−b)² ≥ ¾(a+b)². ∎

Lemma 5 is purely geometric and holds for *any* labelled set S; it is applied below to S = Λ_K and to S = Λ. (Statement (b) for R = 7 is also confirmed by an exhaustive search, `core.geometric_check(7)`, which is called by `tests/equivariance.py`.)

### 4.3 The tiling is well defined (Lemma 6)

σ(S) is the set of children of the sites of a labelled set S; Λ_K := σ^K({0}) (the origin carries the label (0,+1)); lev(y) := min{K : y ∈ Λ_K}.

**Lemma 6.** (a) For every K the 4^K sequences of child indices ("formal sites of depth K") have pairwise distinct positions. Hence Λ_K is a labelled set with 4^K sites, Λ_K ⊂ Λ_{K+1}, and Λ = ⋃_K Λ_K satisfies σ(Λ) = Λ.
 (b) Every y ∈ Λ has exactly one pair (u, j) with u ∈ Λ, j ∈ {0,1,2,3}, y = c_j(u). Write π(y) := u and j(y) := j (for the origin, π(0) = 0 and j(0) = 0).
 (c) lev(π(y)) = lev(y) − 1 for y ≠ 0; thus every y reaches 0 after exactly lev(y) parent steps.

*Proof.* **Step 1 (collisions are local).** Suppose c_j(u) = c_{j′}(u′) for two different pairs (u, j) ≠ (u′, j′). If u = u′ then ε_j = ε_{j′}, which is impossible for j ≠ j′ because the four vectors 0, D_{r+3}, D_{r+2}, D_{r+5} are distinct. So u ≠ u′, and 2(u − u′) = ε_{j′}(u′) − ε_j(u) has hexd ≤ 2; hence hexd(u − u′) = 1. *Only lattice-adjacent sites can claim the same point.*
**Step 2 (finite closure; computation C1a).** Let W⁰ be the radius-7 window of the origin in Λ₈ (`closure.origin_window`; it asserts that Λ₉ gives the same window). Let 𝒲 be the least set of G-classes of windows that contains the class of W⁰ and is closed under the four maps ν_j (`closure.child_windows_plain`; well defined on classes by Lemma 4). The program found |𝒲| = 3542 and inflated every W ∈ 𝒲 without any collision (`core.inflate` asserts `f not in fine` for each pair of children of sites of W).
**Step 3 (induction on K).** A pattern P is *dominated* by a window W (P ⊑ W) if every site of P is a site of W with the same label. Claim S(K): Λ_K has pairwise distinct positions, and for every y ∈ Λ_K there are g ∈ G and W ∈ 𝒲 such that g·P₇(y; Λ_K) ⊑ W. S(0) holds ({0} ⊑ W⁰). Assume S(K).
 – Suppose two formal children of Λ_K coincide. They are children of two adjacent sites u ≠ u′ of Λ_K (Step 1). Take g, W for u: the site u′ lies in P₇(u; Λ_K), so W contains the site g·(u′ − u), and the children of u and u′, transported by g (Lemma 4: g commutes with taking children), are children of two sites of W that coincide. This contradicts Step 2. So Λ_{K+1} = σ(Λ_K) has 4^{K+1} distinct positions.
 – For y′ = c_j(u) ∈ Λ_{K+1}, Lemma 5(b) (with S = Λ_K) gives P₇(y′; Λ_{K+1}) = ν_j(P₇(u; Λ_K)). With g, W for u, g·P₇(y′; Λ_{K+1}) is dominated by ν_j(W), a representative of an element of 𝒲 (Step 2). This is S(K+1).
**Step 4.** Λ₀ ⊂ Λ₁ (0 = c₀(0) with the same label) and σ is monotone for inclusion of labelled sets, so Λ_K ⊂ Λ_{K+1} for all K; with Step 3 this gives (a). Let y ∈ Λ with lev(y) = K ≥ 1. Then y ∈ σ(Λ_{K−1}), so y = c_j(u) with u ∈ Λ_{K−1}; if lev(u) ≤ K − 2 then y ∈ σ(Λ_{K−2}) = Λ_{K−1}, contradicting lev(y) = K. So lev(u) = K − 1. If also y = c_{j′}(u′) with u′ ∈ Λ, choose K′ ≥ K with u′ ∈ Λ_{K′}; then (u, j) and (u′, j′) are two formal children of Λ_{K′} at the same position, so they coincide by (a). This proves (b) and (c). ∎

### 4.4 All windows of the tiling (Lemma 7)

**Lemma 7.** (i) Λ ∩ B₁(0) = {0, D₃, D₂, D₅} and Λ ∩ B₇(0) ⊂ Λ₄. In particular W⁰ (computed in Λ₈) is the radius-7 window of the origin in Λ.
 (ii) For every y ∈ Λ the G-class of W₇(y) belongs to 𝒲.

*Proof.* (i) For y ∈ Λ, y ≠ 0, with parent u = π(y): 2u = y − ε, so hexd(u) ≤ ⌊(hexd(y) + 1)/2⌋; in particular hexd(y) ≤ 1 implies hexd(u) ≤ 1. The children of 0 are 0, D₃, D₂, D₅ (labels (0,+), (0,+), (5,+), (2,−)), so Λ₁ = {0, D₃, D₂, D₅}. The twelve children of D₃ [label (0,+)], D₂ [label (5,+)], D₅ [label (2,−)] sit at (−2,0), (−3,0), (−3,1), (−1,−1); (−2,2), (−3,3), (−2,3), (−2,1); (2,−2), (3,−3), (2,−3), (2,−1) — all with hexd ≥ 2 — so Λ₂ ∩ B₁ = Λ₁ ∩ B₁. If some site of B₁ had level K ≥ 2, its parent chain would stay in B₁ and contain a site of B₁ of level exactly 2 (Lemma 6(c)), contradicting Λ₂ ∩ B₁ = Λ₁ ∩ B₁. Hence Λ ∩ B₁ = {0, D₃, D₂, D₅}, all of level ≤ 1. Now let y ∈ Λ ∩ B₇ with K = lev(y) ≥ 3. Along the parent chain hexd ≤ 7, 4, 2, 1, so the ancestor of level K − 3 lies in B₁, hence K − 3 ≤ 1, i.e. K ≤ 4.
(ii) Induction on lev(y). For y = 0 this is (i). If y = c_j(p) with p = π(y) (lev(p) = lev(y) − 1), then W₇(y) = ν_j(W₇(p)) by Lemma 5(b) (S = Λ), whose class lies in 𝒲 because 𝒲 is closed under the ν_j (and Lemma 4). ∎

### 4.5 Seeds, cells, coverage

*Seeds.* Let 𝒮 ⊂ 𝒲 be the set of classes of windows W₇(z), z ∈ Λ, such that either
 (A) z has no same-species lattice neighbour and all six lattice neighbours of z are sites; or
 (B) z has exactly one same-species lattice neighbour z′, z′ has no same-species lattice neighbour other than z, and all lattice neighbours of z and of z′ are sites.
These conditions involve only P₂(z; Λ) ⊂ W₇(z), so they are decided inside the window. Computation C1b (`closure.seed_windows`, applied to the 3542 classes of 𝒲): **103** classes of type (A) and **22** of type (B). By (H0) and Lemma 7(ii): for chain A the window of the site of C₀ is in 𝒮; for chain B the windows of both sites of C₀ are in 𝒮.

*The map Φ.* For a window W = W₇(z) let Φ(W) be the set of radius-7 windows of the sites c₀(z), c₁(z), c₂(z) (the vertices of T_z) and c₃(w) for every lattice neighbour w of z that is a site and has species different from that of z. By Lemma 5(b) (the triangle vertices) and 5(d) (the c₃(w)), Φ(W) is a function of W — the program asserts the sufficient Euclidean condition |f| + 7 ≤ 7√3 − 1 ≈ 11.12 for the centre f of every output window, which with Lemma 5(e) puts the whole output window inside the determined region B₁₃(2z) — and by Lemma 4 it is G-equivariant, so Φ acts on G-classes (`closure.new_member_windows`).

*Cells.* 𝕎 is the least G-invariant set of windows that contains 𝒮 and is closed under Φ (Computation C1c, `closure.member_closure`): **|𝕎| = 3223** classes. (Admitting also the seeds adjacent to holes of the tiling gives 3440; see §8.) All elements of 𝕎 are windows of Λ, hence in 𝒲.

**Lemma 8 (coverage).** For both chains, every k ≥ 0 and every z ∈ C_k, the G-class of W₇(z) belongs to 𝕎.

*Proof.* Induction on k. For k = 0 this is the statement about seeds. Let y ∈ C_{k+1}. By (H1) and (H2), either y ∈ T_z for some z ∈ C_k, or y = c₃(w) with w a lattice neighbour of some z ∈ C_k; in the second case species(w) = −species(y) = −species(C_{k+1}) = −species(z). In both cases W₇(y) is one of the windows of Φ(W₇(z)), which lies in 𝕎 by the induction hypothesis and closure. ∎

(Φ deliberately returns the windows of *all* fourth children of other-species neighbours, not only those that turn out to belong to the cluster of T_z: whether a given c₃(w) belongs to that cluster can depend on information far outside any fixed window, and nothing below needs to know.)

### 4.6 Local data, types, pieces

Fine positions are measured from 2z (z the centre of a cell): the point 2z + f is called f.

*Rows owned by z.* O(z) := T_z ∪ {c₃(w) : w a lattice neighbour of z that is a site, species(w) ≠ species(z)}. For f ∈ O(z): |f| ≤ 3 (triangle vertices: |f| ≤ 1; c₃(w) = 2w + ε₃, so f = 2(w − z) + ε₃). By (H1), (H2): **every vertex of every cluster C_k (k ≥ 1) of both chains lies in O(z) for some z ∈ C_{k−1}** (for an extra y = c₃(w) with w adjacent to v ∈ C_{k−1}: species(w) = −species(y) = −species(v), so w is of the other species).

*Graph balls.* For a site c ∈ Λ and s ≥ 0, Ball_s(c) is the set of sites reachable from c by at most s steps of 𝒢. If c lies in a cluster C then Ball_s(c) ⊂ C, and Ball_s(c) is read off from the labelled pattern P_s(c; Λ).

**The lower-edge form.** Let k ≥ 1, C = C_k. For sites x and v put β_x(v) := #{g ∈ T_v : g is a same-species lattice neighbour of x}. For x ∈ C and v ∈ C_{k−1}: (Au_v)_x = β_x(v). If β_x(v) ≠ 0 for some v ∈ Λ, then T_v contains a same-species neighbour g of x; g ∈ C is a triangle vertex (child index ≠ 3), so by (H1) and uniqueness of parents (Lemma 6(b)) v = π(g) ∈ C_{k−1}. Hence, for x, y ∈ C,

  (Q₁)_{xy} = q(x,y) := (9/8)δ_{xy} − [x ~ y] + 2 Σ_{v ∈ Λ} β_x(v) β_y(v),

which is nonzero only if dist_𝒢(x,y) ≤ 3 (x ~ g₁, y ~ g₂ with g₁, g₂ ∈ T_v: g₁ = g₂ or g₁ ~ g₂).
*Type and frame.* τ_L(c) := (j(c), [P₂(π(c); Λ)]), the child index together with the G-class of the radius-2 window of the parent. Let (k_c, s_c) ∈ G carry P₂(π(c); Λ) to its class representative (it is unique: every non-identity element of G changes the label of the centre, so windows have trivial stabilisers). The *canonical frame* of c is y ↦ rot^{k_c}(y − c).
*Shape property (S_L).* The set Off(c) := {rot^{k_c}(y − c) : y ∈ Ball_2(c)} depends only on τ_L(c). (By Lemma 5(b) with R = 2, P₂(c; Λ) = ν_{j(c)}(P₂(π(c); Λ)); with Lemma 4 this is a function of τ_L(c) and covariant under G; Ball_2(c) is read off from P₂(c; Λ).)
*Certificate.* A finite set 𝒯_L of types and, for each τ ∈ 𝒯_L, a symmetric positive definite matrix G_τ indexed by Off(τ). For a cluster C: B_c := Ball_2(c), G_c := G_{τ_L(c)} transported to B_c by the canonical frame (this requires τ_L(c) ∈ 𝒯_L), F̂ := Σ_{c ∈ C} ι_cᵀ G_c ι_c, R := Q₁ − F̂.

**The upper-edge form.** For sites v, w ∈ Λ and i ≥ 0 let W_i(v,w) be the number of walks of length i in 𝒢 from a vertex of T_v to a vertex of T_w. For v, w ∈ C_{k−1}, u_vᵀAⁱu_w = W_i(v,w) (walks stay inside the cluster), so with C20 = (−23, −164, −300, 32, 228, 80), the coefficients of 20p (§3.4),

  20·M_{vw} = Σ_{i=0}^{5} C20_i · W_i(v,w) =: q_U(v,w).

Let H_Λ be the graph on Λ in which v ~ w (v ≠ w) iff a vertex of T_v and a vertex of T_w are same-species lattice neighbours.
 (i) v ~_H w implies hexd(v − w) = 1 [f = 2v + ε ∈ T_v and g = 2w + ε′ ∈ T_w adjacent: 2(v − w) = (f − g) − ε + ε′ has hexd ≤ 3, so hexd(v − w) ≤ 3/2, i.e. = 1 as v ≠ w];
 (ii) whether v ~_H w depends only on the labels and the relative position of v and w (the vertices of T_v are determined by the position and label of v);
 (iii) C_{k−1} is closed under H_Λ-adjacency [if v ∈ C_{k−1} and v ~_H w, then a vertex of T_w is a same-species neighbour of a vertex of T_v ⊂ C_k, so T_w ⊂ C_k consists of triangle vertices of C_k, and w ∈ C_{k−1} by (H1) and uniqueness of parents].
So for v ∈ C_{k−1} the H_Λ-ball Ball^H_s(v) lies in C_{k−1} and is read off from P_s(v; Λ).
*Type and frame.* τ_U(c) := [P₄(c; Λ)], the class of the radius-4 window of the coarse site c, with canonical frame y ↦ rot^{k_c}(y − c) (coarse offsets). *Shape property (S_U):* Off^H(c) := {rot^{k_c}(w − c) : w ∈ Ball^H_2(c)} depends only on τ_U(c) (by (i), (ii), Lemma 4).
*Certificate.* Symmetric positive definite G_τ (τ ∈ 𝒯_U) indexed by Off^H(τ); M̂ := Σ_{c ∈ C_{k−1}} ι_cᵀ G_{τ_U(c)} ι_c, R := 20M − M̂.

**Lemma 9 (locality).** Let k ≥ 1, z ∈ C_{k−1}, W = W₇(z), and use fine coordinates f (point 2z + f).
 (a) For |f| ≤ 10: 2z + f ∈ Λ iff f ∈ inflate(W) (the children of the sites of W); label, child index and parent are then those recorded by `inflate(W)`. Consequently, for |f| ≤ 9 the neighbour list `nbrs[f]` of the program is the true list of same-species lattice neighbours of the point.
 (b) (Rows of Q₁.) Let x ∈ O(z). The centres c ∈ Ball_2(x) have |c| ≤ 5 and parents at hexd ≤ 3 from z, so P₂(π(c); Λ) ⊂ W₅(z) ⊂ W; their balls Ball_2(c), and all sites y with q(x,y) ≠ 0 or F̂_{xy} ≠ 0 (graph distance ≤ 4 from x), have |y| ≤ 7; β_y is evaluated from the neighbours of such y (|g| ≤ 8). Everything the program uses for row x is therefore determined by W, via (a).
 (c) (Row of M at z.) Walks of length ≤ 5 from T_z stay within |f| ≤ 6, so the triangles of all endpoints (|f| ≤ 7) lie in the zone and in the program's component `comp` of T_z. The centres c ∈ Ball^H_2(z) have hexd(c − z) ≤ 2 by (i), so P₄(c; Λ) ⊂ B₆(z) ⊂ W; expansions of H-balls happen at sites at H-distance ≤ 3 from z, whose triangle vertices have |f| ≤ 7 and whose neighbours, with their triangles, have |f| ≤ 9 — so all these adjacencies are decided inside the zone. All data of the row are determined by W.

*Proof.* (a) By Lemma 5(e), |f| ≤ 10 ⇒ hexd(f) ≤ 11 ≤ 13, and Lemma 5(c) applies. For |f| ≤ 9 all six lattice neighbours of f have |·| ≤ 10. (b), (c): the bounds follow from |f| ≤ hexd(f), Lemma 5(a) (hexd(π(c) − z) ≤ ⌊(hexd(c) + 1)/2⌋ with hexd(c) ≤ 5) and the facts of §4.6 (i)–(iii). ∎
The program asserts, for every site at which it expands a graph ball, a walk or an H-ball, or evaluates β, that |f| ≤ ZONE − 1 = 9 (`cells.py`, `verify.py`); the bounds above show that this guard is not close to being active (the sites at which anything is expanded or evaluated have |f| ≤ 7).

### 4.7 Lemma 10: the exact check implies Proposition P

**Lemma 10.** Assume the following exact statements.
 (C2) every G_τ (τ ∈ 𝒯_L ∪ 𝒯_U) is symmetric and positive definite.
 (C3L) For every cell W ∈ 𝕎 with centre z and every x ∈ O(z): (i) for every c ∈ Ball_2(x) (computed in the cell), τ_L(c) ∈ 𝒯_L and the canonical frame maps Ball_2(c) onto Off(τ_L(c)); (ii) the row R_{x·} computed from the cell satisfies R_xx − Σ_{y≠x} |R_xy| ≥ ρ_L.
 (C3U) For every cell W ∈ 𝕎 with centre z: (i) for every c ∈ Ball^H_2(z), τ_U(c) ∈ 𝒯_U and the frame maps Ball^H_2(c) onto Off^H(τ_U(c)); (ii) the row of 20M − M̂ at z computed from the cell has dominance margin ≥ ρ_U.
Then, for both chains and every k ≥ 1: Q₁ ⪰ ρ_L·I on ℝ^{C_k} and 20M ⪰ ρ_U·I on ℝ^{C_{k−1}}.

*Proof.* *Lower edge.* Fix k, C = C_k. Let c ∈ C. By §4.6 c ∈ O(z) for some z ∈ C_{k−1}; by Lemma 8, W₇(z) = g·W for some W ∈ 𝕎 and g ∈ G. In the cell W the point x′ := g⁻¹c corresponding to c is an owned row, and the computation of that row examines x′ itself as a centre (x′ ∈ Ball_2(x′)); so (C3L)(i) and Lemma 4 give τ_L(c) ∈ 𝒯_L and the shape property for c. Hence G_c, F̂ and R are defined. Now fix x ∈ C and choose z ∈ C_{k−1}, W, g with x ∈ O(z). By Lemma 9 the row Q₁(x,·) and every ingredient of the row F̂(x,·) — the sets Ball_2(c), the types and canonical frames of the centres c ∈ Ball_2(x), hence the matrices G_{τ_L(c)} and their placement — are obtained from W₇(z) = g·W by explicit rules that commute with g (Lemma 4; the canonical frame is defined through the class representative, so it is G-invariant as a set of offsets). So the row of R at x is the g-image of the row computed in the cell W for g⁻¹x, and its dominance margin is ≥ ρ_L by (C3L)(ii). By (C2) and Lemma 3 (V = C, B_c = Ball_2(c), G_c as above), Q₁ ⪰ ρ_L·I.
*Upper edge.* Identical, with V = C_{k−1}, B_c = Ball^H_2(c) ⊂ C_{k−1} (§4.6 (iii)), Lemma 8 for z ∈ C_{k−1}, and Lemma 9(c). ∎

---

## 5. The exact computations

### 5.1 What is computed

All verdicts are computed with python integers and `fractions.Fraction`. Windows are int8 code vectors (code 0 = no site; code = 1 + 4r + 2·[τ = +1] for a site with label (r,τ); entries in the order of `core.ball(R)`); the 12 elements of G act by permuting the offsets and recoding the labels; the class representative of a window is the lexicographically least of its 12 images (any G-invariant choice would do). Numpy is used only for these int8 arrays.

* **C1 (windows and cells).** (a) `closure.all_sites_closure`: the set 𝒲 of Lemmas 6/7 — **3542** classes. (b) `closure.seed_windows`: **103** classes of type (A), **22** of type (B). (c) `closure.member_closure`: the set 𝕎 — **3223** classes (3440 with the seeds adjacent to holes).
* **C2 (pieces).** `verify.Certificate.check_pieces`: each of the **896** pieces for Q₁ (sizes 1…15, entries integer multiples of 2⁻³⁰) and each of the **1028** pieces for 20M (sizes 1…13, entries integer multiples of 2⁻²⁴) has the shape of its offset list, is symmetric (entrywise equality of integers) and is positive definite by exact Gaussian elimination without pivoting (all pivots > 0; the smallest pivots are 10⁻⁶ resp. 0.0102).
* **C3 (rows).** For every cell W ∈ 𝕎 (`cells.Cell`: one inflation step of W) and every owned row x ∈ O(z) (for M: the row of z), `verify.row_margin_L` / `row_margin_U` compute exactly the row of R, check the shape properties (i) of Lemma 10 and return the dominance margin. The minima over all rows of all cells are
 **ρ_L = 85844911/1073741824 = 0.0799487…** (18473 rows in 3223 cells; 19691 rows in 3440 cells) and
 **ρ_U = 696063391/8388608 = 82.97722…** in units of 20M, i.e. M ⪰ 696063391/167772160 ≈ **4.14886**·I (3223 rows; 3440 with all seeds).

By Lemma 10 this proves Proposition P.

### 5.2 Pseudocode of the verification (what `verify.py` does)

```
for each cell W in 𝕎 (centre z):
    cell := inflate(W)                       # labelled fine sites around 2z with roles (child index, parent offset)
    owned := T_z ∪ { c3(w) : w neighbour of z, species(w) ≠ species(z) }
    for x in owned:                                         # lower edge
        acc := {}
        for c in Ball_2(x):                                 # exactly the centres whose ball contains x
            τ := ( j(c), canonical radius-2 window of the parent of c ); k := rotation to the canonical frame
            check  sorted( rot^k(y − c) : y in Ball_2(c) ) == Off_τ          # shape condition
            for y in Ball_2(c):  acc[y] += G_τ[ idx(x), idx(y) ]
        for y in keys(acc) ∪ Ball_3(x):   R[y] := q(x,y) − acc[y]            # exact Fractions
        margin(x) := R[x] − Σ_{y≠x} |R[y]|
    rho_L(W) := min_x margin(x)
    (upper edge) row of 20M at z by walks of length ≤ 5 from T_z; centres c in the H-ball of radius 2 of z; same bookkeeping with H-balls
assert all pieces symmetric and positive definite, and min_W rho(W) > 0
```

### 5.3 What each check in the code guarantees

| code check | guarantees |
|---|---|
| `core.inflate`: `assert f not in fine` | no collision among the children of the sites of a window (Lemma 6, Step 2) |
| `core.build_patch`: `assert q not in new` | no collision in Λ_K for K ≤ 9 (used for the origin window and for the census) |
| `closure.origin_window`: Λ₈ and Λ₉ give the same window | consistency check of Lemma 7(i) |
| `closure.new_member_windows`: `assert euclid(f) + R <= √3·R − 1` | every output window lies in the determined region of the input window (Lemma 5(c), (e)) |
| `closure.seed_windows` | the seed conditions (A), (B) are decided on the window |
| `cells.Cell.dist_ball`, `Cell.B`, `cells.m_row`, `cells.hball(..., inside)`: `euclid ≤ ZONE − 1` | all neighbour lists used are complete (Lemma 9(a)) |
| `verify.Certificate.type_of`: `assert hexd(o) + Rt <= 7` | the window used for a type lies inside the cell window |
| `verify.row_margin_*`: `assert tuple(…) == p['offs']` | shape condition (i) of Lemma 10; also that the type exists (a missing type raises `KeyError`) |
| `verify.Certificate.check_pieces` | all G_τ symmetric, of the right shape, and exactly positive definite (C2) |
| `verify.main`: `assert rhoL > 0 and rhoU > 0` | row dominance with positive margin on all cells (C3(ii)) |
| `counts.main`: `assert len(seqA) == 1` etc. | all seeds of a chain give the same sizes |

### 5.4 How the pieces were found (irrelevant for correctness)

By a semidefinite program (`tools/gen_constraints.py`, `tools/solve_sdp.py`, Appendix A): unknowns G_τ ⪰ 0 for all types, linear equations Σ_c G_{τ(c)}[x,y] = F_{xy} − ε·[x = y] for all rows of all cells (i.e. R = εI exactly), maximise ε (numerically, with SCS). The numerical solution is then shifted by the smallest multiple of the identity that makes every piece positive definite plus 10⁻⁶ (lower edge) / 10⁻² (upper edge, units of 20M), and rounded to multiples of 2⁻³⁰ / 2⁻²⁴. The rounded pieces leave tiny nonzero off-diagonal entries in R which the row-dominance check absorbs. The optimal slack for the shipped lower-edge certificate was ε* = 0.0803, and the margin certified after rounding is 0.0799; for the upper edge the certified margin is 82.977. Statistics of the shipped certificates are in Appendix B.8, their fingerprints in Appendix C.

---

## 6. Validation (not used in the proof)

All of this uses a code path that does not share the cell machinery: `tests/cluster_numpy.py` builds the actual clusters by brute-force inflation of a finite region (numpy); only the tiling rules are common.

* **Census** (B.6): the finite clusters of the level-9 patch all of whose lattice neighbours are sites have exactly the sizes {1, 2, 4, 10, 16, 42, 63, 165, 245, 941, 3583} — the sizes of the two chains — and no others.
* **Sizes** (B.9): chain A: 1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605 (k ≤ 9); chain B: 2, 10, 42, 165, 632, 2396, 9049, 34135, 128737 (k ≤ 8). They coincide with the substitution of §7.
* **Structure** (B.7): H1 and H2 hold for all these clusters.
* **Eigenvalues** (B.4): dense eigensolver for chain A up to k = 6 (|C₆| = 3583) and chain B up to k = 5 (|C₅| = 2396). In every case #{eigenvalues ≥ 9/8} = |C_{k−1}| and no eigenvalue lies in [9/8, 23/20]. Over all these levels (rounded values from B.4) the largest eigenvalue below the window is ≈ 1.007, the smallest above it is ≈ 1.255 or larger (decreasing slowly with k, consistent with the proven bound > 1.15), λ_min(Q₁) ≈ 0.080 or larger and λ_min(20M) ≈ 102.9 or larger — consistent with the certified bounds 0.0799 and 82.977. (A float is not a proof; this is only a sanity check of the exact argument.)
* **Global exact assembly** (B.4, last two numbers on each line): for these actual clusters the certificate pieces are assembled into F̂ globally (exact Fractions, with the types computed from the actual labelled sites, not from cells) and the Gershgorin margin of R = F − F̂ is recomputed directly on the whole cluster: it is ≥ the margins certified on the cells (0.0799…, 82.977…) in every case — equal to them for the smallest examples found. This checks Lemmas 9/10 and the shape condition on thousands of actual vertices independently of the cells.
* **Equivariance** (B.5): the children map on windows commutes with all 12 group elements (200 random windows × 11 non-identity elements × 4 children), and `geometric_check(7)` confirms Lemma 5(b) exhaustively.

**Hand check, chain A, k = 1.** C₁ has 4 vertices (the three vertices of one triangle and one extra); it is the "paw" graph (degrees 1, 2, 2, 3), with characteristic polynomial λ⁴ − 4λ² − 2λ + 1 = (λ + 1)(λ³ − λ² − 3λ + 1) and eigenvalues 2.1701, 0.3111, −1, −1.4812. For g(λ) = λ³ − λ² − 3λ + 1: g(9/8) = −1135/512 < 0, g(23/20) = −18013/8000 < 0, g(1/2) = −5/8 < 0 < g(0) = 1, g(2) = −1 < 0 < g(3) = 10; so exactly one root (in (2,3)) exceeds 9/8, as predicted (m = |C₀| = 1), and none lies in [9/8, 23/20].

---

## 7. Gap labels (proof of Theorem B)

**Finite k.** By Theorem A, for every E ∈ [λ₁, λ₂] the number of eigenvalues of A(C_k) below E is |C_k| − |C_{k−1}|; hence N_k(E) = 1 − |C_{k−1}|/|C_k|. [P] The values for k = 1,…,7 are in the statement of Theorem B (output B.3).

**The size recurrence.** Let N_k(τ) be the number of sites of C_k whose cell (the G-class of the radius-7 window) is τ ∈ 𝕎. Each site z of C_k contributes the three vertices of T_z to C_{k+1}, and each extra c₃(w) of C_{k+1} is attributed to exactly one site of C_k, namely the lattice neighbour z ∈ C_k of w that minimises (dir(w → z) − r_w) mod 6 (dir = index of the lattice direction, r_w = label index of w). The cell of z determines the cells of all vertices attributed to it, which gives a non-negative integer matrix 𝕄 on the 3223 cells with N_{k+1} = 𝕄N_k (`counts.children_types`, `counts.sizes`). The test "T_z lies in the cluster of the triangle in question" is carried out by breadth-first search inside the determined zone; this is the only place where cluster membership is decided locally (it can only fail by missing a connection that leaves the zone).
For each of the 103 chain-A seeds the Krylov space span{N₀, 𝕄N₀, 𝕄²N₀, …} has **exact dimension 12**, and for each chain-B seed (a pair of matching cells) dimension 13 (`counts.krylov_dimension`, exact Fractions); hence s_k := |C_k| = 1ᵀ𝕄ᵏN₀ satisfies a linear recurrence of order ≤ 12 (resp. 13) for all k. Let t_k be the sequence defined from s₀,…,s₅ by the order-6 recurrence
  s_{k+6} = 8s_{k+5} − 22s_{k+4} + 25s_{k+3} − 7s_{k+2} − 5s_{k+1} + 2s_k.
Then s − t satisfies a recurrence of order ≤ 19 (least common multiple of the two characteristic polynomials), so it vanishes identically as soon as it vanishes for 19 consecutive terms; `counts.main` asserts that all seeds of a chain give the same 61 terms s₀,…,s₆₀ and that these satisfy the order-6 recurrence, i.e. s_k = t_k for k ≤ 60. Hence s = t for all k, for both chains. [X, exact integers; initial values 1, 4, 16, 63, 245, 941 (A) and 2, 10, 42, 165, 632, 2396 (B)]
**Caveat.** That N_{k+1} = 𝕄N_k describes the *true* clusters for all k is the "substitution model". It rests on the local (zone-limited) decision of cluster membership in `counts.children_types`; it reproduces the brute-force sizes of §6 for k ≤ 9 (A) and k ≤ 8 (B) [V], but I have not proved that no extra is ever mis-attributed for larger k. Theorem A does **not** depend on this bookkeeping (the cells use all fourth children of other-species neighbours).

**Characteristic polynomial and irreducibility** (`tools/algebra.py`, output B.10; all exact). The characteristic polynomial of the recurrence is x⁶ − 8x⁵ + 22x⁴ − 25x³ + 7x² + 5x − 2 = (x − 1)·f(x) with f(x) = x⁵ − 7x⁴ + 15x³ − 10x² − 3x + 2. Modulo 3, f ≡ x⁵ − x⁴ − x² − 1; it has no root in 𝔽₃ (values 2, 1, 2 at 0, 1, 2) and is not divisible by any monic polynomial of degree 2 (for instance the remainders modulo the three monic irreducible quadratics x² + 1, x² + x + 2, x² + 2x + 2 are x − 1, 2, x − 1; the program tests all 3 + 9 monic candidates of degree 1 and 2). A quintic with no factor of degree ≤ 2 is irreducible, so f is irreducible over 𝔽₃, hence over ℚ (f is monic).

**Roots.** f(3) = −16 < 0 < f(4) = 22, f(−1/2) = −43/32 < 0 < f(−2/5) = 1408/3125, and f(39/100) > 0 > f(2/5) = −28/3125: f has real roots θ ∈ (3,4), r₁ ∈ (−1/2, −2/5), r₂ ∈ (39/100, 2/5). By Sturm's theorem (exact computation) f has exactly three real roots, so the two remaining roots are non-real, z and z̄. The product of the roots of f is −f(0) = −2, so |z|² = 2/(θ·|r₁|·r₂) < 2/(3·(2/5)·(39/100)) = 500/117 < 4.28 and |z| < 2.07 < 3 < θ. Thus θ is the unique root of f of maximal modulus, it is simple, and every root of (x − 1)f(x) other than θ has modulus < 3 (the root 1 included). All six roots are distinct (f irreducible, f(1) = −2 ≠ 0).

**Limit.** Hence s_k = cθᵏ + O(2.07ᵏ) (the root 1 contributes a constant). Since every site of C_{k−1} contributes a triangle, s_k ≥ 3·s_{k−1} and s_k ≥ 3ᵏ s₀; this excludes c = 0, so s_{k−1}/s_k → 1/θ and N_k(E) → ℓ := 1 − 1/θ.
Numerically (`tools/algebra.py`, exact bisection): θ = 3.7725380988570505404160270370041892583877, 1/θ = 0.2650735324059326665550111967750136168646, ℓ = 1 − 1/θ = 0.7349264675940673334449888032249863831354. The number 1/θ is a root of 2y⁵ − 3y⁴ − 10y³ + 15y² − 7y + 1, and ℓ is a root of 2ℓ⁵ − 7ℓ⁴ − 2ℓ³ + 13ℓ² − 9ℓ + 2 (verified symbolically in `tools/algebra.py`). The observed ratios s_{k−1}/s_k (0.26494 for A at k = 9, 0.26515 for B at k = 8) are close to 1/θ.

---

## 8. Generality: chains that start next to holes

Singleton clusters that have a lattice neighbour which is not a site of the tiling (a hole) are finite clusters too, but they are "on the boundary" and are not part of the statement. There are 61 window classes of such seeds (164 − 103); they generate three further families of chains, with sizes (substitution model, output B.11)
 (1, 4, 16, 62, 236, 889, 3336, 12515, …) [9 classes], (1, 4, 15, 56, 210, 789, 2966, 11152, …) [36 classes], (1, 4, 15, 55, 201, 737, 2719, 10093, …) [16 classes].
Adding the seeds to the closure enlarges 𝕎 to 3440 cells; the *same* certificates verify on all 3440 cells (`python -m ww.verify --all-seeds`, output B.2). Hence Proposition P — and, whenever H1/H2 hold for such a chain, Theorem A with the same λ₁, λ₂ — also holds for them. (H1/H2 were not checked for these chains by brute force.)

---

## 9. Assessment: assumptions, risks, what is not proved

**Assumptions [H].** The proof assumes H0–H2, i.e. the description of the chains in the problem statement. H1 and H2 were verified by brute force for k ≤ 9 (A) and k ≤ 8 (B), but are not re-proved for all k here; I took them as part of the definition of the chains, as the problem does. The identification of "the two chains" with the chains that start from a single site or from an adjacent pair of equal-species sites with all neighbours present is my reading, supported by the census of §1.3 (no other finite cluster with all neighbours present occurs in Λ₉). If the intended chains were different clusters, Theorem A as proved here would not cover them (but see §8).

**Computer-assisted steps [X].** C1–C3 cannot be checked by hand. Their statements are explicit (§4.7, §5); the code is in Appendix A (about 650 lines for the part that matters: `core`, `closure`, `cells`, `verify`). The verdict depends on (i) the code implementing the definitions above and (ii) Python's exact integer/rational arithmetic. I mitigated (i) by assertions that encode the lemmas (§5.3), by an independent brute-force implementation (§6) and by recomputing the certified inequality globally on real clusters. I did not formally verify the code.

**Structure of the argument that a judge may want to probe.**
1. Lemmas 1–3 (§3): short and elementary; Lemma 1 holds for P positive semidefinite (a remark gives an invertible variant).
2. Reduction to Proposition P (§3.5): the only inputs are n₊(P) = m and dim S = m, both consequences of H1 (disjoint triangles).
3. Lemma 5 (geometry), Lemma 6 (the tiling is collision free: the one slightly delicate induction, which uses the closure computation to exclude two adjacent sites claiming the same edge midpoint at any level), Lemma 7.
4. Coverage (Lemma 8), locality (Lemma 9) and the certificate lemma (Lemma 10): bookkeeping with explicit distance bounds, plus the symmetry argument that makes a finite number of cells sufficient.
5. The numbers: ρ_L ≈ 0.0799 and ρ_U/20 ≈ 4.149 are small compared with the entries of the forms, as the hints warned; that is why a sum-of-squares certificate checked in exact arithmetic is used rather than a diagonal-dominance argument on the forms themselves (which fails badly: on C₄ of chain A the minimal Gershgorin margins of Q₁ and of 20M are −112.9 and −3.2·10⁵).

**Not claimed.** Optimality of the window (numerically the gap is wider: the largest eigenvalue below it is ≈ 1.007 and the smallest above it ≳ 1.25 for k ≤ 6); that every finite cluster of the tiling lies on these chains beyond the census of §1.3; the substitution model for all k (the limit in Theorem B depends on it; Theorem A does not).

---

## Appendix A — source code (verbatim)

All files are Python (3.13 was used) with numpy; `tools/solve_sdp.py` additionally needs scipy, cvxpy and scs; `ww/counts.py` uses sympy for one optional line (irreducibility). Run commands from the package root (the directory containing `ww/`). Paths refer to the repository directory `walton-whittaker-gap/`.

### A.1 `ww/core.py` — tiling rules, windows, group action

```python
"""Tiling rules, label windows and the 12-element symmetry group of the charged Walton-Whittaker tiling.

Lattice sites are pairs (a, b) = a*D0 + b*D1 with D_j = exp(i*pi*j/3).  A site carries a label (r, tau),
r in Z/6, tau = +-1.  Inflation of a site at p with label (r, tau) gives four children
    2p,            2p + D_{r+3}   labelled (r, tau)          (children 0, 1)
    2p + D_{r+2}                  labelled (r-1, tau)        (child 2)      -> children 0,1,2 form a triangle
    2p + D_{r+5}                  labelled (r+2, -tau)       (child 3, the "flipped" one).
"""
import math
import numpy as np

D = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]


def hexd(a, b):
    """lattice (hexagonal) distance of the offset (a, b)"""
    return max(abs(a), abs(b), abs(a + b))


def rot(p):
    """rotation by +60 degrees: D_j -> D_{j+1}"""
    return (-p[1], p[0] + p[1])


def rotk(p, k):
    for _ in range(k):
        p = rot(p)
    return p


def ball(R):
    return [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if hexd(a, b) <= R]


def euclid(p):
    a, b = p
    return math.sqrt((a + 0.5 * b) ** 2 + 0.75 * b * b)


def child_pos(o, r, j):
    """fine position of child j of the coarse site at offset o with label index r"""
    if j == 0:
        c = (0, 0)
    elif j == 1:
        c = D[(r + 3) % 6]
    elif j == 2:
        c = D[(r + 2) % 6]
    else:
        c = D[(r + 5) % 6]
    return (2 * o[0] + c[0], 2 * o[1] + c[1])


def child_label(r, t, j):
    return ((r % 6, t), (r % 6, t), ((r + 5) % 6, t), ((r + 2) % 6, -t))[j]


def build_patch(K, label=(0, 1)):
    """K-fold inflation of a single site at the origin: dict position -> label"""
    s = {(0, 0): label}
    for _ in range(K):
        new = {}
        for p, (r, t) in s.items():
            for j in range(4):
                q = child_pos(p, r, j)
                assert q not in new, "collision"
                new[q] = child_label(r, t, j)
        s = new
    return s


# ---------------------------------------------------------------------------------------------
# Windows: label content of the hexagonal ball of radius R around a site, as an int8 code vector.
# code 0 = no site; code = 1 + 4 r + 2 [tau = +1] + m  (m is a spare bit, always 0 here).
# ---------------------------------------------------------------------------------------------
def enc(r, t, m=0):
    return 1 + 4 * r + 2 * (1 if t == 1 else 0) + m


def dec(c):
    c -= 1
    return c >> 2, (1 if (c >> 1) & 1 else -1), c & 1


class Frame:
    """windows of radius R and the action of the group G = Z_6 (rotations) x Z_2 (species flip)"""

    def __init__(self, R):
        self.R = R
        self.offs = ball(R)
        self.idx = {o: i for i, o in enumerate(self.offs)}
        inv_rot = lambda p: (p[0] + p[1], -p[0])
        self.perm = np.array([self.idx[inv_rot(o)] for o in self.offs], dtype=np.int64)
        self.rotcode = np.zeros(25, dtype=np.int8)
        self.flipcode = np.zeros(25, dtype=np.int8)
        for c in range(1, 25):
            r, t, m = dec(c)
            self.rotcode[c] = enc((r + 1) % 6, t, m)
            self.flipcode[c] = enc(r, -t, m)

    def rotate(self, v):
        return self.rotcode[v[self.perm]]

    def flip(self, v):
        return self.flipcode[v]

    def canon(self, v):
        """(key, k, s): lexicographically least of the 12 images; k = number of 60-degree rotations used"""
        best = None
        cur = v
        for k in range(6):
            for s, a in ((0, cur), (1, self.flip(cur))):
                b = a.tobytes()
                if best is None or b < best[0]:
                    best = (b, k, s)
            cur = self.rotate(cur)
        return best

    def to_dict(self, v):
        return {o: dec(int(v[i])) for i, o in enumerate(self.offs) if v[i]}

    def from_dict(self, d):
        v = np.zeros(len(self.offs), dtype=np.int8)
        for o, (r, t, m) in d.items():
            if o in self.idx:
                v[self.idx[o]] = enc(r, t, m)
        return v

    def window_of(self, lab, p):
        """window of a labelled point set (dict position -> (r, t)) around position p"""
        v = np.zeros(len(self.offs), dtype=np.int8)
        for i, o in enumerate(self.offs):
            l = lab.get((p[0] + o[0], p[1] + o[1]))
            if l is not None:
                v[i] = enc(l[0], l[1])
        return v


def inflate(d):
    """d: dict coarse offset -> (r, t, m).  Returns dict fine position -> (r, t, j, parent offset).
    Asserts that no two children coincide (a collision would need two adjacent sites to claim the midpoint of their
    common edge); applied to every window of the closure this proves, by induction on the level, that the tiling is
    collision free: a collision at level K+1 would be seen inside the window of one of the two parents."""
    fine = {}
    for o, (r, t, m) in d.items():
        for j in range(4):
            lab = child_label(r, t, j)
            f = child_pos(o, r, j)
            assert f not in fine, "collision of two children"
            fine[f] = (lab[0], lab[1], j, o)
    return fine


def geometric_check(R):
    """Window of radius R of a child is determined by the radius-R window of the parent:
    every parent of a site within distance R of a child lies within distance R of the parent."""
    cj = [(0, 0)] + list(D)
    for c0 in cj:
        for o in ball(R + 3):
            if hexd(*o) <= R:
                continue
            for c in cj:
                if hexd(2 * o[0] + c[0] - c0[0], 2 * o[1] + c[1] - c0[1]) <= R:
                    return False
    return True
```

### A.2 `ww/closure.py` — all windows, seeds, cells

```python
"""Finite set of inflation cells.

1. all-sites closure: every radius-R window of the tiling, obtained from the window of the origin by iterating
   the children map (a window of a child is a function of the window of its parent, Lemma 5 of FULL_PROOF.md);
2. seeds: windows centred at the member of a singleton cluster (chain A) or of a two-site cluster (chain B);
3. member closure: windows of the sites of C_k that are vertices of the triangle of a site z of C_{k-1} or flipped fourth
   children of neighbours of z, as a function of the window of z.  The least set containing the seeds and closed under
   this map is the set of cells (up to the symmetry group G).
"""
import math
import numpy as np
from .core import D, Frame, child_pos, euclid, inflate, build_patch

ZONE = 10.0   # fine sites farther than this (euclidean) from the inflated centre are not determined by a radius-7 window


# ---------------------------------------------------------------------------------------------
# 1. all windows of the tiling
# ---------------------------------------------------------------------------------------------
def origin_window(F, K=8, label=(0, 1)):
    patch = build_patch(K, label)
    w = {o: (patch[o][0], patch[o][1], 0) for o in F.offs if o in patch}
    patch2 = build_patch(K + 1, label)
    assert w == {o: (patch2[o][0], patch2[o][1], 0) for o in F.offs if o in patch2}, "origin window not stable"
    return w


def child_windows_plain(F, d):
    """windows of the four children of the centre of the window d"""
    fine = inflate(d)
    r0 = d[(0, 0)][0]
    out = []
    for j in range(4):
        cj = child_pos((0, 0), r0, j)
        w = {}
        for f, v in fine.items():
            o = (f[0] - cj[0], f[1] - cj[1])
            if o in F.idx:
                w[o] = (v[0], v[1], 0)
        out.append(w)
    return out


def all_sites_closure(F):
    d0 = origin_window(F)
    k0, _, _ = F.canon(F.from_dict(d0))
    seen = {k0: np.frombuffer(k0, dtype=np.int8).copy()}
    frontier = [seen[k0]]
    while frontier:
        nf = []
        for v in frontier:
            for w in child_windows_plain(F, F.to_dict(v)):
                k, _, _ = F.canon(F.from_dict(w))
                if k not in seen:
                    seen[k] = np.frombuffer(k, dtype=np.int8).copy()
                    nf.append(seen[k])
        frontier = nf
    return list(seen.values())


# ---------------------------------------------------------------------------------------------
# 2. seeds
# ---------------------------------------------------------------------------------------------
def same_species_nbrs(d, o):
    t = d[o][1]
    return [(o[0] + dl[0], o[1] + dl[1]) for dl in D
            if (o[0] + dl[0], o[1] + dl[1]) in d and d[(o[0] + dl[0], o[1] + dl[1])][1] == t]


def seed_windows(F, plain, interior_only=True):
    """('A', v): window centred at a singleton cluster; ('B', v): window centred at a site of a two-site cluster.
    interior_only: every lattice neighbour of every site of the seed cluster is a site of the tiling
    ("off the boundary"); the other seeds start chains next to holes of the tiling."""
    seeds = []
    z = (0, 0)
    for v in plain:
        d = F.to_dict(v)
        nb = same_species_nbrs(d, z)
        if len(nb) == 0:
            if interior_only and not all(dl in d for dl in D):
                continue
            seeds.append(('A', v))
        elif len(nb) == 1 and same_species_nbrs(d, nb[0]) == [z]:
            if interior_only and not all((s[0] + dl[0], s[1] + dl[1]) in d for s in (z, nb[0]) for dl in D):
                continue
            seeds.append(('B', v))
    return seeds


# ---------------------------------------------------------------------------------------------
# 3. member closure
# ---------------------------------------------------------------------------------------------
def component_from(fine, starts, zone=ZONE):
    """same-species connected component (nearest-neighbour graph) of the fine sites 'starts'"""
    comp = set(starts)
    fr = list(starts)
    while fr:
        nf = []
        for f in fr:
            sp = fine[f][1]
            for dl in D:
                g = (f[0] + dl[0], f[1] + dl[1])
                v = fine.get(g)
                if v is None or g in comp or v[1] != sp or euclid(g) > zone:
                    continue
                comp.add(g)
                nf.append(g)
        fr = nf
    return comp


def new_member_windows(F, d, bfs=False):
    """Windows of the possible new members associated with the centre z of d: the three triangle children of z and the
    flipped fourth children w_3 of the other-species neighbours w of z.  Every vertex of the cluster C_{k+1} that is
    not a vertex of a triangle T_z is such an extra w_3 with w adjacent to some z of C_k; the cell does not need to know
    whether w_3 really belongs to the cluster of z (that may be decided far outside the window), so all of them are
    returned (bfs=True returns only those found in the cluster of T_z inside the determined zone; used for counting)."""
    fine = inflate(d)
    z = (0, 0)
    r0, t0, _ = d[z]
    cands = [child_pos(z, r0, j) for j in range(3)]
    for q in D:
        if q in d and d[q][1] != t0:
            cands.append(child_pos(q, d[q][0], 3))
    if bfs:
        comp = component_from(fine, cands[:3])
        cands = [f for f in cands if f in comp]
    outs = []
    for f in cands:
        assert euclid(f) + F.R <= math.sqrt(3) * F.R - 1.0 + 1e-9, "output window not determined by input window"
        w = {}
        for o in F.offs:
            g = (f[0] + o[0], f[1] + o[1])
            v = fine.get(g)
            if v is not None:
                w[o] = (v[0], v[1], 0)
        outs.append(w)
    return outs


def member_closure(F, seeds, bfs=False):
    seen = {}
    frontier = []
    for _, v in seeds:
        k, _, _ = F.canon(v)
        if k not in seen:
            seen[k] = np.frombuffer(k, dtype=np.int8).copy()
            frontier.append(seen[k])
    while frontier:
        nf = []
        for v in frontier:
            for w in new_member_windows(F, F.to_dict(v), bfs):
                k, _, _ = F.canon(F.from_dict(w))
                if k not in seen:
                    seen[k] = np.frombuffer(k, dtype=np.int8).copy()
                    nf.append(seen[k])
        frontier = nf
    return list(seen.values())


def cells(R=7, interior_only=True):
    F = Frame(R)
    plain = all_sites_closure(F)
    seeds = seed_windows(F, plain, interior_only)
    return F, plain, seeds, member_closure(F, seeds)
```

### A.3 `ww/cells.py` — inflation cell, the two local forms

```python
"""An inflation cell: one inflation step applied to the radius-7 window of a site z of C_{k-1}.

It gives the labelled fine sites around 2z (inside the region determined by the window), their roles (child index, parent),
the triangles (children 0,1,2 of one parent), the extras (child 3) and the same-species nearest-neighbour graph.  The
two quadratic forms below are computed from it exactly.
"""
from fractions import Fraction as Fr
from .core import D, child_pos, inflate, euclid
from .closure import component_from, ZONE

# 20 * p5(x) with p5(x) = (2x^2+4x+1)^2 (x - 23/20) = 4x^5 + 57/5 x^4 + 8/5 x^3 - 15 x^2 - 41/5 x - 23/20
C20 = [-23, -164, -300, 32, 228, 80]          # coefficients of x^0 .. x^5


class Cell:
    def __init__(self, d):
        self.d = d
        self.fine = inflate(d)
        z = (0, 0)
        r0, t0, _ = d[z]
        self.tri = [child_pos(z, r0, j) for j in range(3)]
        self.comp = component_from(self.fine, self.tri)      # cluster of the triangle of z (inside the determined zone)
        # nearest-neighbour graph of ALL fine sites inside the zone: the graph ball around a site of a cluster is the
        # ball of this graph, whether or not the site is known (inside the zone) to belong to the cluster of T_z
        zone = [f for f in self.fine if euclid(f) <= ZONE]
        zs = set(zone)
        self.nbrs = {f: [(f[0] + dl[0], f[1] + dl[1]) for dl in D
                         if (f[0] + dl[0], f[1] + dl[1]) in zs and self.fine[(f[0] + dl[0], f[1] + dl[1])][1] == self.fine[f][1]]
                     for f in zone}
        self.tris = {}                      # parent offset -> fine sites of its triangle (cluster of T_z only)
        for f in self.comp:
            v = self.fine[f]
            if v[2] != 3:
                self.tris.setdefault(v[3], []).append(f)
        # fine sites whose rows of Q1 are checked in this cell: the triangle of z and EVERY flipped fourth child of an
        # other-species neighbour of z (whether it belongs to the cluster of z or to another one)
        self.owned = list(self.tri)
        for q in D:
            if q in d and d[q][1] != t0:
                self.owned.append(child_pos(q, d[q][0], 3))
        self._B = {}

    def dist_ball(self, x, rad):
        dist = {x: 0}
        fr = [x]
        for s in range(rad):
            nf = []
            for y in fr:
                assert euclid(y) <= ZONE - 1, "graph ball leaves the determined zone"
                for g in self.nbrs[y]:
                    if g not in dist:
                        dist[g] = s + 1
                        nf.append(g)
            fr = nf
        return dist

    def B(self, x):
        """B_x[v] = (A u_v)_x = number of neighbours of x inside the triangle of the coarse site v"""
        if x not in self._B:
            assert euclid(x) <= ZONE - 1, "site too close to the edge of the determined zone"
            out = {}
            for g in self.nbrs[x]:
                v = self.fine[g]
                if v[2] != 3:
                    out[v[3]] = out.get(v[3], 0) + 1
            self._B[x] = out
        return self._B[x]


# ------------------------------------------------------------------------------------------------
# lower edge:  Q1 = 9/8 I - A + 2 sum_v (A u_v)(A u_v)^T   (P = 2 sum_v u_v u_v^T)
# ------------------------------------------------------------------------------------------------
def q1_entry(cell, x, y, pw=2):
    val = Fr(0)
    if x == y:
        val += Fr(9, 8)
    if y in cell.nbrs[x]:
        val -= 1
    Bx, By = cell.B(x), cell.B(y)
    return val + pw * sum(c * By[v] for v, c in Bx.items() if v in By)


# ------------------------------------------------------------------------------------------------
# upper edge:  M = U^T p5(A) U  (U = triangle indicator vectors), rows computed as 20*M
# ------------------------------------------------------------------------------------------------
def coarse_data(cell):
    tri = {po: sites for po, sites in cell.tris.items() if len(sites) == 3}
    site2tri = {f: po for po, sites in tri.items() for f in sites}
    H = {po: set() for po in tri}
    for po, sites in tri.items():
        for f in sites:
            for g in cell.nbrs[f]:
                if g in site2tri and site2tri[g] != po:
                    H[po].add(site2tri[g])
    return tri, site2tri, H


def hball(H, x, rad, inside=None):
    """ball of radius rad in the triangle-adjacency graph H; inside(y) asserts that y is far enough from the zone edge"""
    dist = {x: 0}
    fr = [x]
    for s in range(rad):
        nf = []
        for y in fr:
            assert inside is None or inside(y), "H-ball leaves the determined zone"
            for g in H[y]:
                if g not in dist:
                    dist[g] = s + 1
                    nf.append(g)
        fr = nf
    return dist


def m_row(cell, tri, site2tri, po):
    """20 * M[po, w] for all coarse sites w: walks of length <= 5 starting in the triangle po"""
    a = [{f: 1 for f in tri[po]}]
    for _ in range(5):
        nxt = {}
        for f, v in a[-1].items():
            assert euclid(f) <= ZONE - 1, "walk leaves the determined zone"
            for g in cell.nbrs[f]:
                nxt[g] = nxt.get(g, 0) + v
        a.append(nxt)
    row = {}
    for i, ai in enumerate(a):
        for f, v in ai.items():
            w = site2tri.get(f)
            if w is not None:
                row[w] = row.get(w, 0) + C20[i] * v
    return row
```

### A.4 `ww/verify.py` — exact verification (the program whose output is in B.1)

```python
"""Exact-arithmetic verification of the two local certificates over all inflation cells.

For a quadratic form F on the vertices of a cluster (F = Q1 on fine vertices, F = 20*M on coarse sites) the certificate is
a family of positive definite matrices G_tau (tau = local type of a centre c) such that

        F = sum_c  iota_c^T G_{tau(c)} iota_c  +  R,        R_xx >= sum_{y != x} |R_xy| + rho   for every row x,

where iota_c embeds the ball W_c around the centre c.  Since every term of the sum is positive semidefinite and R is
diagonally dominant (Gershgorin), F >= rho * I.  The row inequality only involves the cell of x, so it is checked
for the finitely many cells.  Everything below is computed with python integers / Fractions.
"""
import argparse
import collections
import gzip
import json
import os
import time
from fractions import Fraction as Fr

import numpy as np

from .core import Frame, enc, euclid, hexd, rotk
from .cells import Cell, q1_entry, coarse_data, hball, m_row
from .closure import cells as build_cells, ZONE

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')


def is_pd(M):
    """exact Gaussian elimination without pivoting: (True, smallest pivot) iff the symmetric matrix M is positive definite"""
    n = len(M)
    A = [row[:] for row in M]
    minp = None
    for k in range(n):
        p = A[k][k]
        if p <= 0:
            return False, p
        minp = p if minp is None or p < minp else minp
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = A[i][k] / p
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return True, minp


class Certificate:
    def __init__(self, path):
        with gzip.open(path, 'rt') as fh:
            raw = json.load(fh)
        self.form, self.Rt, self.bits = raw['form'], raw['Rt'], raw['bits']
        self.FR = Frame(self.Rt)
        self.pieces = {}
        for t in raw['types']:
            key = (t['j'], bytes.fromhex(t['key']))
            self.pieces[key] = dict(offs=tuple(tuple(o) for o in t['offs']), N=t['N'])

    def check_pieces(self):
        """every G_tau symmetric and positive definite (exact); returns the smallest pivot met"""
        S = Fr(1 << self.bits)
        minpiv = None
        for key, p in self.pieces.items():
            n = len(p['N'])
            assert len(p['offs']) == n and all(len(row) == n for row in p['N']), "piece shape"
            assert all(p['N'][i][j] == p['N'][j][i] for i in range(n) for j in range(i)), "piece not symmetric"
            ok, mp = is_pd([[Fr(v) / S for v in row] for row in p['N']])
            assert ok, "piece not positive definite"
            minpiv = mp if minpiv is None or mp < minpiv else minpiv
        return minpiv

    def type_of(self, d, o, j):
        """canonical type of the centre with child index j whose parent sits at offset o of the window d"""
        FR, Rt = self.FR, self.Rt
        assert hexd(*o) + Rt <= 7
        v = np.zeros(len(FR.offs), dtype=np.int8)
        for i, q in enumerate(FR.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        key, k, s = FR.canon(v)
        return (j, key), k

    def piece(self, d, o, j):
        """the piece matrix (integers over 2^bits, canonical frame) of the centre, and the rotation k to that frame"""
        tkey, k = self.type_of(d, o, j)
        return self.pieces[tkey], k


def row_margin_L(cert, cell, x, pw=2):
    S = Fr(1 << cert.bits)
    cache = {}

    def centre(c):
        if c not in cache:
            v = cell.fine[c]
            p, k = cert.piece(cell.d, v[3], v[2])
            ball_ = cell.dist_ball(c, 2)
            coffs = {y: rotk((y[0] - c[0], y[1] - c[1]), k) for y in ball_}
            order = sorted(ball_, key=lambda y: coffs[y])
            assert tuple(coffs[y] for y in order) == p['offs'], "piece window mismatch"
            cache[c] = (p, {y: i for i, y in enumerate(order)})
        return cache[c]

    acc = collections.defaultdict(Fr)
    for c in cell.dist_ball(x, 2):
        p, idx = centre(c)
        ix = idx[x]
        for y, iy in idx.items():
            acc[y] += Fr(p['N'][ix][iy]) / S
    cols = set(acc) | set(cell.dist_ball(x, 3))
    R = {y: q1_entry(cell, x, y, pw) - acc.get(y, Fr(0)) for y in cols}
    return R[x] - sum(abs(v) for y, v in R.items() if y != x)


def row_margin_U(cert, cell, Wrad=2):
    S = Fr(1 << cert.bits)
    tri, site2tri, H = coarse_data(cell)
    z = (0, 0)
    cache = {}
    inside = lambda po: all(euclid(f) <= ZONE - 1 for f in tri[po])

    def centre(c):
        if c not in cache:
            p, k = cert.piece(cell.d, c, -1)
            ball_ = hball(H, c, Wrad, inside)
            coffs = {w: rotk((w[0] - c[0], w[1] - c[1]), k) for w in ball_}
            order = sorted(ball_, key=lambda w: coffs[w])
            assert tuple(coffs[w] for w in order) == p['offs'], "piece window mismatch"
            cache[c] = (p, {w: i for i, w in enumerate(order)})
        return cache[c]

    row = m_row(cell, tri, site2tri, z)
    acc = collections.defaultdict(Fr)
    for c in hball(H, z, Wrad, inside):
        p, idx = centre(c)
        ix = idx[z]
        for w, iw in idx.items():
            acc[w] += Fr(p['N'][ix][iw]) / S
    cols = set(acc) | set(row)
    R = {w: Fr(row.get(w, 0)) - acc.get(w, Fr(0)) for w in cols}
    return R[z] - sum(abs(v) for w, v in R.items() if w != z)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--all-seeds', action='store_true', help='also include the chains that start next to holes of the tiling')
    ap.add_argument('--data', default=DATA)
    args = ap.parse_args(argv)
    t0 = time.time()
    F, plain, seeds, wins = build_cells(7, interior_only=not args.all_seeds)
    nA = sum(1 for c, _ in seeds if c == 'A')
    print(f"all windows of the tiling (radius 7, mod G): {len(plain)}")
    print(f"seeds: {nA} singleton windows (chain A), {len(seeds) - nA} pair windows (chain B)")
    print(f"inflation cells (closure under inflation, mod G): {len(wins)}   [{time.time() - t0:.0f}s]")
    L = Certificate(os.path.join(args.data, 'cert_L.json.gz'))
    U = Certificate(os.path.join(args.data, 'cert_U.json.gz'))
    pL, pU = L.check_pieces(), U.check_pieces()
    print(f"lower edge: {len(L.pieces)} pieces, all symmetric and positive definite (exact), smallest pivot {float(pL):.3g}")
    print(f"upper edge: {len(U.pieces)} pieces, all symmetric and positive definite (exact), smallest pivot {float(pU):.3g}")
    rhoL = rhoU = None
    nrows = 0
    for n, v in enumerate(wins):
        cell = Cell(F.to_dict(v))
        for x in cell.owned:
            r = row_margin_L(L, cell, x)
            rhoL = r if rhoL is None or r < rhoL else rhoL
            nrows += 1
        r = row_margin_U(U, cell)
        rhoU = r if rhoU is None or r < rhoU else rhoU
        if n % 500 == 0:
            print(f"  cell {n}/{len(wins)}  [{time.time() - t0:.0f}s]", flush=True)
    print(f"lower edge: {nrows} rows checked; min row margin rho_L = {rhoL} = {float(rhoL):.6f} > 0")
    print(f"upper edge: {len(wins)} rows checked; min row margin rho_U = {rhoU} = {float(rhoU):.6f} > 0 (units of 20*M)")
    assert rhoL > 0 and rhoU > 0
    print(f"=> Q1 = 9/8 - A + 2 sum_v (A u_v)(A u_v)^T  >=  {float(rhoL):.4f} * I")
    print(f"=> M = U^T p5(A) U  >=  {float(rhoU) / 20:.4f} * I     (u_v = indicator of the triangle of v, not normalised)")
    print(f"VERIFIED in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
```

### A.5 `ww/counts.py` — substitution, sizes, recurrence, gap labels

```python
"""Cluster sizes |C_k| from the substitution on inflation cells, their linear recurrence, and the gap labels.

Every site of C_{k-1} contributes the 3 vertices of its triangle; every extra (flipped fourth child of a neighbouring site
w of the other species) is attributed to the neighbour z* of w, inside the same cluster, with the smallest
(dir(w -> z') - r_w) mod 6.  The cell of z determines the cells of all attributed new members, so the sizes obey a
linear recurrence; the Perron root gives lim |C_{k-1}|/|C_k|, i.e. the integrated density of states above the gap.
"""
import collections
import sys
from fractions import Fraction as Fr

import numpy as np

from .core import D, Frame, child_pos, inflate, enc
from .closure import cells as build_cells, same_species_nbrs, component_from

REC = [8, -22, 25, -7, -5, 2]       # s_{k+6} = 8 s_{k+5} - 22 s_{k+4} + 25 s_{k+3} - 7 s_{k+2} - 5 s_{k+1} + 2 s_k
CHARPOLY = "x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x-1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2)"


def children_types(F, d, type_index):
    """types (indices into the cell list) of the new members attributed to the centre of the cell d"""
    fine = inflate(d)
    z = (0, 0)
    r0, t0, _ = d[z]
    tri = [child_pos(z, r0, j) for j in range(3)]
    comp = component_from(fine, tri)
    cands = list(tri)
    for q in D:
        if q in d and d[q][1] != t0:
            f = child_pos(q, d[q][0], 3)
            if f not in comp:
                continue
            rw = d[q][0]
            members = []
            for dl in D:
                zz = (q[0] + dl[0], q[1] + dl[1])
                if zz in d and d[zz][1] == t0 and all(child_pos(zz, d[zz][0], j) in comp for j in range(3)):
                    members.append(zz)
            assert z in members
            zstar = min(members, key=lambda zz: (D.index((zz[0] - q[0], zz[1] - q[1])) - rw) % 6)
            if zstar == z:
                cands.append(f)
    out = []
    for f in cands:
        w = {}
        for o in F.offs:
            g = (f[0] + o[0], f[1] + o[1])
            if g in fine:
                w[o] = (fine[g][0], fine[g][1], 0)
        k, _, _ = F.canon(F.from_dict(w))
        out.append(type_index[k])
    return out


def sizes(trans, start, K):
    vec = dict(start)
    s = [sum(vec.values())]
    for _ in range(K):
        nv = collections.Counter()
        for t, c in vec.items():
            for ch in trans[t]:
                nv[ch] += c
        vec = nv
        s.append(sum(vec.values()))
    return s


def krylov_dimension(trans, start, maxdim=60):
    """exact dimension of span{N_0, m N_0, m^2 N_0, ...}: the vectors N_k of cell counts live in this space, so every
    linear functional of N_k (e.g. |C_k|) obeys a linear recurrence of at most this order"""
    vec = {t: Fr(c) for t, c in start.items()}
    basis = []
    for _ in range(maxdim):
        v = dict(vec)
        for piv, b in basis:
            c = v.get(piv, 0)
            if c:
                for k, x in b.items():
                    v[k] = v.get(k, 0) - c * x
                    if v[k] == 0:
                        del v[k]
        if not v:
            return len(basis)
        piv = min(v)
        c = v[piv]
        basis.append((piv, {k: x / c for k, x in v.items()}))
        nv = {}
        for t, c in vec.items():
            for ch in trans[t]:
                nv[ch] = nv.get(ch, 0) + c
        vec = nv
    return None


def main():
    F, plain, seeds, wins = build_cells(7, interior_only=True)
    type_index = {v.tobytes(): i for i, v in enumerate(wins)}
    trans = [children_types(F, F.to_dict(v), type_index) for v in wins]
    A_seeds = [v for c, v in seeds if c == 'A']
    seqA = {tuple(sizes(trans, {type_index[v.tobytes()]: 1}, 60)) for v in A_seeds}
    assert len(seqA) == 1, "all singleton seeds give the same sizes"
    sA = list(seqA.pop())
    kdA = {krylov_dimension(trans, {type_index[v.tobytes()]: 1}, maxdim=200) for v in A_seeds}
    print(f"exact Krylov dimension of the substitution started from a chain-A seed (all {len(A_seeds)} seeds): {sorted(kdA)}  "
          f"(=> |C_k| obeys a linear recurrence of order <= {max(kdA)} for ALL k; the minimal one has order 6)")
    # chain B: a pair seed = two member windows that see each other (matched through their radius-6 windows)
    F6 = Frame(6)

    def sub6(d, o):
        v = np.zeros(len(F6.offs), dtype=np.int8)
        for i, q in enumerate(F6.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        return F6.canon(v)[0]

    B = [v for c, v in seeds if c == 'B']
    info = []
    for v in B:
        d = F.to_dict(v)
        info.append((v, sub6(d, (0, 0)), sub6(d, same_species_nbrs(d, (0, 0))[0])))
    seqB, kdB = set(), set()
    for i, (v, c6, p6) in enumerate(info):
        for j, (w, c6j, p6j) in enumerate(info):
            if c6j == p6 and p6j == c6:
                start = collections.Counter([type_index[v.tobytes()], type_index[w.tobytes()]])
                seqB.add(tuple(sizes(trans, start, 60)))
                kdB.add(krylov_dimension(trans, start, maxdim=200))
    assert len(seqB) == 1
    sB = list(seqB.pop())
    print(f"exact Krylov dimension of the substitution started from a chain-B seed (all matched pairs): {sorted(kdB)}")
    assert max(kdA | kdB) + 6 <= len(sA), "too few terms to conclude that the recurrence holds for all k"
    for name, s in (('A', sA), ('B', sB)):
        assert all(s[k + 6] == sum(c * s[k + 5 - i] for i, c in enumerate(REC)) for k in range(len(s) - 6))
        print(f"chain {name}: |C_k| = {s[:10]} ...   (satisfies the recurrence with characteristic polynomial {CHARPOLY})")
    quintic = [1, -7, 15, -10, -3, 2]
    theta = max(r.real for r in np.roots(quintic) if abs(r.imag) < 1e-12)
    print(f"Perron root theta = {theta:.12f},  lim |C_(k-1)|/|C_k| = 1/theta = {1 / theta:.12f},  1 - 1/theta = {1 - 1 / theta:.12f}")
    print("finite-k labels 1 - |C_(k-1)|/|C_k|  (A):", [str(1 - Fr(sA[k - 1], sA[k])) for k in range(1, 8)])
    print("finite-k labels 1 - |C_(k-1)|/|C_k|  (B):", [str(1 - Fr(sB[k - 1], sB[k])) for k in range(1, 8)])
    try:
        import sympy as sp
        x = sp.symbols('x')
        print("quintic irreducible over Q:", sp.Poly(sum(c * x ** (5 - i) for i, c in enumerate(quintic)), x).is_irreducible)
    except ImportError:
        pass


if __name__ == '__main__':
    main()
```

### A.6 `tools/algebra.py` — exact algebra of the gap labels

```python
"""Exact algebra behind the gap labels (Theorem B of FULL_PROOF.md); python integers / Fractions / Decimal only.

  * the characteristic polynomial of the size recurrence factors as (x - 1) f(x), f(x) = x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2;
  * f is irreducible over Q: no root and no quadratic factor modulo 3 (exhaustive test of all monic polynomials of degree 1, 2);
  * Sturm's theorem: f has exactly three real roots; they are isolated by sign changes at rational points;
  * hence every root other than theta (the root in (3, 4)) has modulus < 3 < theta;
  * theta, 1/theta, 1 - 1/theta to 40 digits (exact bisection) and the polynomials they satisfy;
  * the recurrence reproduces the brute-force sizes of tests/sizes.py.
usage: python -m tools.algebra
"""
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
from itertools import product

REC = [8, -22, 25, -7, -5, 2]                      # s_{k+6} = 8 s_{k+5} - 22 s_{k+4} + 25 s_{k+3} - 7 s_{k+2} - 5 s_{k+1} + 2 s_k
F = [2, -3, -10, 15, -7, 1]                        # f, coefficients of x^0 .. x^5


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pev(p, x):
    r = 0
    for c in reversed(p):
        r = r * x + c
    return r


def trim(p):
    while p and p[-1] == 0:
        p = p[:-1]
    return p


def pmod(a, b, mod=None):
    a = [Fr(c) for c in a] if mod is None else [c % mod for c in a]
    b = trim([Fr(c) for c in b] if mod is None else [c % mod for c in b])
    a = trim(a)
    while len(a) >= len(b):
        if mod is None:
            f = a[-1] / b[-1]
        else:
            f = a[-1] * pow(b[-1], -1, mod) % mod
        sh = len(a) - len(b)
        for i, c in enumerate(b):
            a[sh + i] = a[sh + i] - f * c if mod is None else (a[sh + i] - f * c) % mod
        a = trim(a)
    return a


def sturm_count(p):
    """number of distinct real roots of p (Sturm's theorem, exact)"""
    seq = [[Fr(c) for c in p], [Fr(i * c) for i, c in enumerate(p)][1:]]
    while True:
        r = pmod(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(signs):
        s = [x for x in signs if x != 0]
        return sum(1 for a, b in zip(s, s[1:]) if a != b)

    at_minus = [(1 if len(q) % 2 == 1 else -1) * (1 if q[-1] > 0 else -1) for q in seq]    # sign at -infinity
    at_plus = [1 if q[-1] > 0 else -1 for q in seq]                                       # sign at +infinity
    return changes(at_minus) - changes(at_plus)


def main():
    getcontext().prec = 60
    P6 = [-REC[5], -REC[4], -REC[3], -REC[2], -REC[1], -REC[0], 1]      # x^6 - 8x^5 + ... - 2, ascending coefficients
    assert P6 == pmul([-1, 1], F), "char poly = (x - 1) f(x)"
    print("x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x - 1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2): verified")

    # irreducibility mod 3
    f3 = [c % 3 for c in F]
    print("f mod 3 =", f3, "(ascending coefficients)")
    assert all(pev(f3, a) % 3 != 0 for a in range(3)), "root mod 3"
    for deg in (1, 2):
        for low in product(range(3), repeat=deg):
            g = list(low) + [1]
            assert trim(pmod(f3, g, 3)), f"f has the factor {g} modulo 3"
    print("f mod 3 has no monic factor of degree 1 or 2 (all 3 + 9 candidates tested): f is irreducible over F_3, hence over Q")

    # real roots
    n_real = sturm_count(F)
    assert n_real == 3
    print("Sturm: f has exactly", n_real, "real roots")
    pts = [Fr(-1, 2), Fr(-2, 5), Fr(39, 100), Fr(2, 5), Fr(3), Fr(4)]
    vals = [pev(F, x) for x in pts]
    print("f at -1/2, -2/5, 39/100, 2/5, 3, 4:", [str(v) for v in vals])
    assert vals[0] < 0 < vals[1] and vals[2] > 0 > vals[3] and vals[4] < 0 < vals[5]
    # product of the roots = -f(0) = -2: the two non-real roots z, conj(z) satisfy |z|^2 = 2 / (theta |r1| r2)
    bound = Fr(2) / (3 * Fr(2, 5) * Fr(39, 100))
    print(f"non-real roots: |z|^2 = 2/(theta |r1| r2) < 2/(3 * 2/5 * 39/100) = {bound} = {float(bound):.4f} < 9, so |z| < 3 < theta")
    assert bound < 9

    # theta by exact bisection
    lo, hi = Fr(3), Fr(4)
    for _ in range(140):
        mid = (lo + hi) / 2
        if pev(F, mid) < 0:
            lo = mid
        else:
            hi = mid
    th = (lo + hi) / 2

    def dec(q):
        return Decimal(q.numerator) / Decimal(q.denominator)

    print("theta        =", f"{dec(th):.40f}")
    print("1/theta      =", f"{dec(1 / th):.40f}")
    print("1 - 1/theta  =", f"{dec(1 - 1 / th):.40f}")

    # polynomials satisfied by y = 1/theta and l = 1 - 1/theta
    y_poly = list(reversed(F))                     # y^5 f(1/y)
    print("y = 1/theta is a root of", " + ".join(f"({c})y^{i}" for i, c in enumerate(y_poly)), "= 0")
    # l = 1 - y  <=>  y = 1 - l : substitute into y_poly (Horner with polynomial arithmetic in l)
    res = [0]
    for c in reversed(y_poly):
        res = pmul(res, [1, -1])
        res[0] += c
    res = trim(res)
    print("l = 1 - 1/theta is a root of", " + ".join(f"({c})l^{i}" for i, c in enumerate(res)), "= 0")
    assert abs(float(pev(res, 1 - 1 / th))) < 1e-12 and abs(float(pev(y_poly, 1 / th))) < 1e-12

    # the recurrence against the brute-force sizes of tests/sizes.py
    A = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605]
    B = [2, 10, 42, 165, 632, 2396, 9049, 34135, 128737]
    for name, s in (('A', A), ('B', B)):
        pred = list(s[:6])
        for k in range(6, len(s)):
            pred.append(sum(c * pred[k - 1 - i] for i, c in enumerate(REC)))
        assert pred == s, name
        print(f"chain {name}: the recurrence started from s_0..s_5 reproduces the brute-force sizes s_6..s_{len(s) - 1}")
        print(f"   ratios s_(k-1)/s_k:", [f"{float(Fr(s[k - 1], s[k])):.5f}" for k in range(1, len(s))])


if __name__ == '__main__':
    main()
```

### A.7 `tools/gen_constraints.py` — constraints of the piece SDP

```python
"""Generate the linear constraints of the piece SDP over all inflation cells.

Lower edge (L):  pieces G_tau sit on the ball of graph radius 2 around a fine vertex c; type tau = (child index j of c,
                 canonical radius-Rt window of the parent of c).
Upper edge (U):  pieces sit on the ball of radius 2 (triangle adjacency graph H) around a coarse site c; type tau =
                 canonical radius-Rt window of c.
For every row x of every cell the equation  sum_c G_{tau(c)}[x, y] = F_xy  (minus a slack eps on the diagonal) is one
constraint; constraints with the same left-hand side and right-hand side are merged.

usage: python -m tools.gen_constraints L 2 out.pkl [--all-seeds]       (Rt = 2 for L, 4 for U)
"""
import argparse
import collections
import pickle
import sys
from fractions import Fraction as Fr

import numpy as np

from ww.core import Frame, enc, hexd, rotk
from ww.cells import Cell, q1_entry, coarse_data, hball, m_row
from ww.closure import cells as build_cells


class TypeBook:
    def __init__(self, Rt):
        self.Rt, self.FR = Rt, Frame(Rt)
        self.ids, self.keys, self.windows = {}, [], []

    def type_of(self, d, o, j):
        FR = self.FR
        assert hexd(*o) + self.Rt <= 7
        v = np.zeros(len(FR.offs), dtype=np.int8)
        for i, q in enumerate(FR.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        key, k, s = FR.canon(v)
        full = (j, key)
        if full not in self.ids:
            self.ids[full] = len(self.ids)
            self.keys.append(full)
            self.windows.append(None)
        return self.ids[full], k

    def set_window(self, tid, win):
        if self.windows[tid] is None:
            self.windows[tid] = win
        assert self.windows[tid] == win, "the type does not determine the piece window: increase Rt"


def rows_L(tb, cell):
    cinfo = {}

    def centre(c):
        if c not in cinfo:
            v = cell.fine[c]
            tid, k = tb.type_of(cell.d, v[3], v[2])
            ball_ = cell.dist_ball(c, 2)
            coffs = {y: rotk((y[0] - c[0], y[1] - c[1]), k) for y in ball_}
            order = sorted(ball_, key=lambda y: coffs[y])
            tb.set_window(tid, tuple(coffs[y] for y in order))
            cinfo[c] = (tid, {y: i for i, y in enumerate(order)})
        return cinfo[c]

    for x in cell.owned:
        contrib = collections.defaultdict(list)
        for c in cell.dist_ball(x, 2):
            tid, idx = centre(c)
            for y, iy in idx.items():
                contrib[y].append((tid, idx[x], iy))
        for y in set(contrib) | set(cell.dist_ball(x, 3)):
            yield q1_entry(cell, x, y), tuple(sorted(contrib.get(y, []))), y == x


def rows_U(tb, cell, Wrad=2):
    tri, site2tri, H = coarse_data(cell)
    z = (0, 0)
    cinfo = {}

    def centre(c):
        if c not in cinfo:
            tid, k = tb.type_of(cell.d, c, -1)
            ball_ = hball(H, c, Wrad)
            coffs = {w: rotk((w[0] - c[0], w[1] - c[1]), k) for w in ball_}
            order = sorted(ball_, key=lambda w: coffs[w])
            tb.set_window(tid, tuple(coffs[w] for w in order))
            cinfo[c] = (tid, {w: i for i, w in enumerate(order)})
        return cinfo[c]

    row = m_row(cell, tri, site2tri, z)
    contrib = collections.defaultdict(list)
    for c in hball(H, z, Wrad):
        tid, idx = centre(c)
        for w, iw in idx.items():
            contrib[w].append((tid, idx[z], iw))
    for w in set(contrib) | set(row):
        assert not (row.get(w, 0) != 0 and w not in contrib), "uncovered nonzero entry"
        yield Fr(row.get(w, 0)), tuple(sorted(contrib.get(w, []))), w == z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('form', choices=['L', 'U'])
    ap.add_argument('Rt', type=int)
    ap.add_argument('out')
    ap.add_argument('--all-seeds', action='store_true')
    a = ap.parse_args()
    F, plain, seeds, wins = build_cells(7, interior_only=not a.all_seeds)
    tb = TypeBook(a.Rt)
    sigs = {}
    for n, v in enumerate(wins):
        cell = Cell(F.to_dict(v))
        for rhs, lst, isdiag in (rows_L(tb, cell) if a.form == 'L' else rows_U(tb, cell)):
            sigs[(lst, rhs, isdiag)] = 1
    print(f"{a.form}: {len(wins)} cells, {len(tb.ids)} types, {len(sigs)} constraints")
    pickle.dump(dict(Rt=a.Rt, keys=tb.keys, windows=tb.windows, sigs=list(sigs)), open(a.out, 'wb'))


if __name__ == '__main__':
    main()
```

### A.8 `tools/solve_sdp.py` — SDP solve, shift, rounding, export

```python
"""Solve the piece SDP (numerically, SCS) and turn the solution into a rational certificate.

    maximise eps   s.t.   sum_c G_{tau(c)}[x,y] = F_xy - eps*[x=y]   (all constraints),   G_tau positive semidefinite.

The numerical solution is only a candidate: every piece is shifted by the smallest multiple of the identity that makes
it positive definite and rounded to a multiple of 2^-bits; ww.verify then checks everything exactly.

usage: python -m tools.solve_sdp constraints.pkl out.json.gz FORM [--tol 1e-4] [--scale 1] [--bits 30] [--margin 1e-6]
"""
import argparse
import gzip
import json
import pickle

import numpy as np
import scipy.sparse as sp
import cvxpy as cp


def solve(windows, sigs, tol, scale):
    sizes = [len(w) for w in windows]
    offs = [0]
    for n in sizes:
        offs.append(offs[-1] + n * n)
    merged = {}
    for lst, rhs, isdiag in sigs:
        nl = tuple(sorted((t, min(i, j), max(i, j)) for t, i, j in lst))
        merged[(nl, isdiag)] = rhs
    rows, cols, b, diag = [], [], [], []
    for r, ((nl, isdiag), rhs) in enumerate(merged.items()):
        for t, i, j in nl:
            rows.append(r)
            cols.append(offs[t] + j * sizes[t] + i)
        b.append(float(rhs) / scale)
        diag.append(1.0 if isdiag else 0.0)
    Aeq = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(merged), offs[-1]))
    G = [cp.Variable((n, n), PSD=True) for n in sizes]
    eps = cp.Variable()
    vec = cp.hstack([cp.vec(g, order='F') for g in G])
    prob = cp.Problem(cp.Maximize(eps), [Aeq @ vec + eps * np.array(diag) == np.array(b)])
    prob.solve(solver='SCS', eps=tol, max_iters=200000)
    return prob.status, float(eps.value) * scale, [g.value * scale for g in G]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('constraints')
    ap.add_argument('out')
    ap.add_argument('form', choices=['L', 'U'])
    ap.add_argument('--tol', type=float, default=1e-4)
    ap.add_argument('--scale', type=float, default=1.0)
    ap.add_argument('--bits', type=int, default=30)
    ap.add_argument('--margin', type=float, default=1e-6)
    a = ap.parse_args()
    data = pickle.load(open(a.constraints, 'rb'))
    status, eps, pieces = solve(data['windows'], data['sigs'], a.tol, a.scale)
    print(f"SDP status {status}, optimal slack eps* = {eps:.6g}")
    S = 1 << a.bits
    types = []
    for (j, key), win, G in zip(data['keys'], data['windows'], pieces):
        G = (G + G.T) / 2
        G = G + (max(0.0, -np.linalg.eigvalsh(G).min()) + a.margin) * np.eye(len(G))
        N = [[0] * len(G) for _ in G]
        for i in range(len(G)):
            for k in range(i, len(G)):
                N[i][k] = N[k][i] = int(round(float(G[i, k]) * S))
        types.append(dict(j=j, key=key.hex(), offs=[list(o) for o in win], N=N))
    with gzip.open(a.out, 'wt') as fh:
        json.dump(dict(form=a.form, Rt=data['Rt'], bits=a.bits, types=types), fh, separators=(',', ':'))
    print(f"{len(types)} pieces written to {a.out}")


if __name__ == '__main__':
    main()
```

### A.9 `tools/cert_stats.py` — statistics and fingerprints of the certificate files

```python
"""Summary statistics of the shipped certificates (sizes of pieces, magnitude of entries) and fingerprints of the data."""
import collections
import gzip
import hashlib
import json
import os

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')
for form in ('L', 'U'):
    path = os.path.join(DATA, f'cert_{form}.json.gz')
    raw = json.load(gzip.open(path, 'rt'))
    types = raw['types']
    sizes = collections.Counter(len(t['N']) for t in types)
    mx = max(abs(v) for t in types for row in t['N'] for v in row)
    print(f"certificate {form}: type radius Rt = {raw['Rt']}, entries are integers over 2^{raw['bits']}, {len(types)} pieces")
    print(f"  piece sizes (size: number of pieces): {dict(sorted(sizes.items()))}")
    if form == 'L':
        print(f"  pieces by child index j of the centre: {dict(sorted(collections.Counter(t['j'] for t in types).items()))}")
    print(f"  largest |entry| = {mx / 2 ** raw['bits']:.6g}")
    print(f"  sha256 {hashlib.sha256(open(path, 'rb').read()).hexdigest()}  (gzip file, {os.path.getsize(path)} bytes)")
    small = min(types, key=lambda t: (len(t['N']), t['j']))
    big = max(types, key=lambda t: len(t['N']))
    for name, t in (('smallest', small), ('largest', big)):
        S = 2 ** raw['bits']
        print(f"  example {name} piece: j = {t['j']}, {len(t['N'])}x{len(t['N'])}, offsets {t['offs'][:4]}{'...' if len(t['offs']) > 4 else ''}")
        for row in t['N'][:3]:
            print("     ", [round(v / S, 5) for v in row[:6]], '...' if len(row) > 6 else '')
```

### A.10 `tests/cluster_numpy.py` — independent brute-force construction of the actual clusters

```python
"""Independent construction of the actual clusters C_k by brute-force inflation of a finite labelled region (numpy).

Used only by the cross-checks: it does not share any code with the cell machinery of ww/ except the tiling rules.
A region = all tiling sites within lattice distance RAD of the current cluster; inflating every site of the region and
taking the same-species connected component of the new triangles gives the next cluster and its region
(RAD >= 6 is self-consistent, see the comment in inflate_cluster).
"""
import numpy as np

from ww.core import D, build_patch

DX = np.array([d[0] for d in D])
DY = np.array([d[1] for d in D])
OFF = 1 << 30


def key(a, b):
    return (a.astype(np.int64) + OFF) * (1 << 32) + (b.astype(np.int64) + OFF)


def unkey(k):
    return (k >> 32) - OFF, (k & ((1 << 32) - 1)) - OFF


def children_arr(a, b, r, t):
    n = len(a)
    idx = np.arange(n)
    out = []
    for j in range(4):
        if j == 0:
            da = np.zeros(n, dtype=np.int64); db = np.zeros(n, dtype=np.int64); rr, tt = r, t
        elif j == 1:
            da = DX[(r + 3) % 6]; db = DY[(r + 3) % 6]; rr, tt = r, t
        elif j == 2:
            da = DX[(r + 2) % 6]; db = DY[(r + 2) % 6]; rr, tt = (r + 5) % 6, t
        else:
            da = DX[(r + 5) % 6]; db = DY[(r + 5) % 6]; rr, tt = (r + 2) % 6, -t
        out.append((2 * a + da, 2 * b + db, rr, tt, idx, np.full(n, j)))
    return [np.concatenate([o[i] for o in out]) for i in range(6)]


def dilate(a, b, steps):
    cur = np.unique(key(a, b))
    allk, frontier = cur, cur
    for _ in range(steps):
        fa, fb = unkey(frontier)
        nk = np.unique(np.concatenate([key(fa + DX[i], fb + DY[i]) for i in range(6)]))
        new = np.setdiff1d(nk, allk, assume_unique=True)
        allk = np.union1d(allk, new)
        frontier = new
    return allk


def component(keys, species, startkeys):
    n = len(keys)
    inc = np.zeros(n, dtype=bool)
    frontier = np.searchsorted(keys, startkeys)
    inc[frontier] = True
    a, b = unkey(keys)
    while len(frontier):
        fa, fb, fs = a[frontier], b[frontier], species[frontier]
        cand = []
        for i in range(6):
            nk = key(fa + DX[i], fb + DY[i])
            p = np.minimum(np.searchsorted(keys, nk), n - 1)
            ok = (keys[p] == nk) & (species[p] == fs) & (~inc[p])
            cand.append(p[ok])
        c = np.unique(np.concatenate(cand)) if cand else np.array([], dtype=np.int64)
        inc[c] = True
        frontier = c
    return inc


def inflate_cluster(a, b, r, t, mask, rad):
    """one inflation step of a region (arrays sorted by key) with the cluster marked by mask.
    Sites within lattice distance 6 of C_k have parents within distance 6 of C_{k-1} (1 + (6+2)/sqrt(3) < 6),
    hence the radius-6 neighbourhood of the new cluster is again generated by the old region."""
    ca, cb, cr, ct, par, ch = children_arr(a, b, r, t)
    keys = key(ca, cb)
    assert len(np.unique(keys)) == len(keys)
    o = np.argsort(keys)
    keys, ca, cb, cr, ct, par, ch = keys[o], ca[o], cb[o], cr[o], ct[o], par[o], ch[o]
    start = keys[mask[par] & (ch != 3)]
    comp = component(keys, ct, start)
    near = dilate(ca[comp], cb[comp], rad)
    inreg = np.isin(keys, near, assume_unique=True)
    return ca[inreg], cb[inreg], cr[inreg], ct[inreg], comp[inreg], par[inreg], ch[inreg]


def seed_region(chain, rad=10, Kpatch=9):
    """a singleton (A) / two-site (B) cluster all of whose lattice neighbours are tiling sites, deep inside a patch,
    together with all tiling sites within distance rad of it"""
    patch = build_patch(Kpatch)
    lim = 2 ** Kpatch // 4

    def same(p):
        t = patch[p][1]
        return [q for q in ((p[0] + d[0], p[1] + d[1]) for d in D) if q in patch and patch[q][1] == t]

    for p in sorted(patch, key=lambda p: -(abs(p[0]) + abs(p[1]))):
        if abs(p[0]) > lim or abs(p[1]) > lim:
            continue
        nb = same(p)
        if chain == 'A' and len(nb) == 0:
            sites = [p]
        elif chain == 'B' and len(nb) == 1 and same(nb[0]) == [p]:
            sites = [p, nb[0]]
        else:
            continue
        if not all((s[0] + d[0], s[1] + d[1]) in patch for s in sites for d in D):
            continue
        pts = np.array(sites)
        near = dilate(pts[:, 0], pts[:, 1], rad)
        na, nb_ = unkey(near)
        sel = [(int(x), int(y)) for x, y in zip(na, nb_) if (int(x), int(y)) in patch]
        a = np.array([q[0] for q in sel]); b = np.array([q[1] for q in sel])
        r = np.array([patch[q][0] for q in sel]); t = np.array([patch[q][1] for q in sel])
        o = np.argsort(key(a, b))
        a, b, r, t = a[o], b[o], r[o], t[o]
        ss = set(sites)
        mask = np.array([(int(x), int(y)) in ss for x, y in zip(a, b)])
        return a, b, r, t, mask
    raise RuntimeError("no seed found")


def chain_levels(chain, L, rad=10):
    """yields (level, region_before, cluster arrays) for k = 1..L"""
    a, b, r, t, mask = seed_region(chain, rad)
    prev = (a, b, r, t, mask)
    for k in range(1, L + 1):
        a, b, r, t, mask, par, ch = inflate_cluster(*prev, rad)
        yield k, prev, (a, b, r, t, mask, par, ch)
        prev = (a, b, r, t, mask)
```

### A.11 `tests/crosscheck.py` — eigenvalue counts, global exact residuals, structure H1/H2

```python
"""Independent numerical / exact cross-checks on the actual clusters C_k (brute-force numpy construction).

For chain A (levels 1..LA) and chain B (levels 1..LB) it
  * builds the cluster graph by brute-force inflation (tests/cluster_numpy.py, no cell machinery),
  * counts the eigenvalues of A(C_k) above 9/8 and in [9/8, 23/20] (dense float eigensolver; a float is NOT the proof),
  * assembles the forms Q1 and 20*M globally and checks their smallest eigenvalues against the certified margins,
  * assembles the certificate pieces globally (exact Fractions) and recomputes the Gershgorin margin of the residual
    R = F - sum_c iota_c^T G_c iota_c on the whole cluster: it must reproduce (be >=) the margins certified on the cells.
  * checks the structural hypotheses H1, H2 of FULL_PROOF.md on the actual clusters: the triangle vertices of C_k are exactly
    the triangles of the sites of C_{k-1}, and every other vertex is the flipped fourth child of a lattice neighbour of a site of C_{k-1}.
Run:  python -m tests.crosscheck [LA LB]              (from the package root)
      python -m tests.crosscheck LA LB structure      (only the structural checks, much faster, deeper levels)
"""
import collections
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

from ww.core import D, child_pos, enc, rotk, Frame
from ww.cells import C20
from ww.verify import Certificate, DATA
from tests.cluster_numpy import chain_levels

RHO_L = Fr(85844911, 1073741824)        # certified by ww.verify for the shipped certificates
RHO_U = Fr(696063391, 8388608)


def cluster_data(prev, cur):
    pa, pb, pr, pt, pmask = prev
    a, b, r, t, mask, par, ch = cur
    cs = [(int(x), int(y)) for x, y, mk in zip(a, b, mask) if mk]
    idx = {p: i for i, p in enumerate(cs)}
    parent = {(int(x), int(y)): int(p) for x, y, mk, p in zip(a, b, mask, par) if mk}
    child = {(int(x), int(y)): int(c) for x, y, mk, c in zip(a, b, mask, ch) if mk}
    plab = {(int(x), int(y)): (int(rr), int(tt)) for x, y, rr, tt in zip(pa, pb, pr, pt)}
    ppos = {i: (int(x), int(y)) for i, (x, y) in enumerate(zip(pa, pb))}
    n_old = int(pmask.sum())
    # structural hypotheses H1, H2 (FULL_PROOF.md, section 1.3)
    old = {(int(x), int(y)) for x, y, mk in zip(pa, pb, pmask) if mk}
    for p in cs:
        w = ppos[parent[p]]
        if child[p] != 3:
            assert w in old, "H1 violated: triangle vertex whose parent is not in C_(k-1)"
        else:
            assert any((w[0] + d[0], w[1] + d[1]) in old for d in D), "H2 violated: extra whose parent is not adjacent to C_(k-1)"
    for v in old:
        assert all(child_pos(v, plab[v][0], j) in idx for j in range(3)), "H1 violated: triangle of a site of C_(k-1) not in C_k"
    adj = [[idx[(p[0] + dl[0], p[1] + dl[1])] for dl in D if (p[0] + dl[0], p[1] + dl[1]) in idx] for p in cs]
    tri = collections.defaultdict(list)
    for p in cs:
        if child[p] != 3:
            tri[parent[p]].append(idx[p])
    return cs, adj, tri, parent, child, plab, ppos, n_old


def dense_forms(adj, tri):
    n = len(adj)
    A = np.zeros((n, n))
    for i, nb in enumerate(adj):
        A[i, nb] = 1
    tids = sorted(tri)
    U = np.zeros((n, len(tids)))
    for k, t in enumerate(tids):
        U[tri[t], k] = 1
    AU = A @ U
    Q1 = 9 / 8 * np.eye(n) - A + 2 * AU @ AU.T
    P = sum(c * np.linalg.matrix_power(A, i) for i, c in enumerate(C20))
    M20 = U.T @ P @ U
    return A, Q1, np.rint(M20), tids


def canon_plain(FR, plab, p):
    v = np.zeros(len(FR.offs), dtype=np.int8)
    for i, o in enumerate(FR.offs):
        l = plab.get((p[0] + o[0], p[1] + o[1]))
        if l is not None:
            v[i] = enc(l[0], l[1], 0)
    return FR.canon(v)


def ball_dist(nbrs, x, rad):
    dist = {x: 0}
    fr = [x]
    for s in range(rad):
        nf = []
        for y in fr:
            for g in nbrs[y]:
                if g not in dist:
                    dist[g] = s + 1
                    nf.append(g)
        fr = nf
    return dist


def margin_L(cert, cs, adj, tri, parent, child, plab, ppos):
    n = len(cs)
    Q = [collections.defaultdict(Fr) for _ in range(n)]
    for i in range(n):
        Q[i][i] += Fr(9, 8)
        for j in adj[i]:
            Q[i][j] -= 1
    for v, members in tri.items():
        Au = collections.Counter()
        for t in members:
            for j in adj[t]:
                Au[j] += 1
        for i, ci in Au.items():
            for j, cj in Au.items():
                Q[i][j] += 2 * ci * cj
    S = Fr(1 << cert.bits)
    Fh = [collections.defaultdict(Fr) for _ in range(n)]
    for c, p in enumerate(cs):
        key, k, s = canon_plain(cert.FR, plab, ppos[parent[p]])
        piece = cert.pieces[(child[p], key)]
        dist = ball_dist(adj, c, 2)
        coffs = {y: rotk((cs[y][0] - p[0], cs[y][1] - p[1]), k) for y in dist}
        order = sorted(dist, key=lambda y: coffs[y])
        assert tuple(coffs[y] for y in order) == piece['offs']
        for a_, y in enumerate(order):
            for b_, w in enumerate(order):
                Fh[y][w] += Fr(piece['N'][a_][b_]) / S
    worst = None
    for i in range(n):
        R = {w: Q[i].get(w, Fr(0)) - Fh[i].get(w, Fr(0)) for w in set(Q[i]) | set(Fh[i])}
        m = R[i] - sum(abs(v) for w, v in R.items() if w != i)
        worst = m if worst is None or m < worst else worst
    return worst


def margin_U(cert, cs, adj, tri, plab, ppos, M20, tids):
    m = len(tids)
    site2tri = {s: k for k, t in enumerate(tids) for s in tri[t]}
    H = [set() for _ in range(m)]
    for k, t in enumerate(tids):
        for s in tri[t]:
            for g in adj[s]:
                if g in site2tri and site2tri[g] != k:
                    H[k].add(site2tri[g])
    S = Fr(1 << cert.bits)
    Mh = [collections.defaultdict(Fr) for _ in range(m)]
    for k, t in enumerate(tids):
        p = ppos[t]
        key, rk, s = canon_plain(cert.FR, plab, p)
        piece = cert.pieces[(-1, key)]
        dist = ball_dist(H, k, 2)
        coffs = {y: rotk((ppos[tids[y]][0] - p[0], ppos[tids[y]][1] - p[1]), rk) for y in dist}
        order = sorted(dist, key=lambda y: coffs[y])
        assert tuple(coffs[y] for y in order) == piece['offs']
        for a_, y in enumerate(order):
            for b_, w in enumerate(order):
                Mh[y][w] += Fr(piece['N'][a_][b_]) / S
    worst = None
    for i in range(m):
        nz = np.nonzero(M20[i])[0]
        R = {w: Fr(int(M20[i, w])) - Mh[i].get(w, Fr(0)) for w in set(nz.tolist()) | set(Mh[i])}
        mm = R[i] - sum(abs(v) for w, v in R.items() if w != i)
        worst = mm if worst is None or mm < worst else worst
    return worst


def main():
    LA = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    LB = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    if len(sys.argv) > 3 and sys.argv[3] == 'structure':
        for chain, Lmax in (('A', LA), ('B', LB)):
            for k, prev, cur in chain_levels(chain, Lmax):
                cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
                assert len(tri) == n_old
                print(f"chain {chain} k={k}: |C_k|={len(cs)}, |C_(k-1)|={n_old}: H1, H2 hold", flush=True)
        return
    L = Certificate(os.path.join(DATA, 'cert_L.json.gz'))
    U = Certificate(os.path.join(DATA, 'cert_U.json.gz'))
    lam1, lam2 = 9 / 8, 23 / 20
    for chain, Lmax in (('A', LA), ('B', LB)):
        sizes = []
        for k, prev, cur in chain_levels(chain, Lmax):
            t0 = time.time()
            cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
            n, m = len(cs), len(tri)
            assert m == n_old, "triangles of C_k <-> sites of C_{k-1}"
            A, Q1, M20, tids = dense_forms(adj, tri)
            w = np.linalg.eigvalsh(A)
            above = int((w >= lam1).sum())
            inwin = int(((w >= lam1) & (w <= lam2)).sum())
            top = w[w > lam2].min()
            bulk = w[w < lam1].max()
            eq1 = np.linalg.eigvalsh(Q1).min()
            em = np.linalg.eigvalsh((M20 + M20.T) / 2).min()
            ml = margin_L(L, cs, adj, tri, parent, child, plab, ppos)
            mu = margin_U(U, cs, adj, tri, plab, ppos, M20, tids)
            ok = (above == m and inwin == 0 and ml >= RHO_L and mu >= RHO_U and eq1 >= float(RHO_L) and em >= float(RHO_U))
            print(f"chain {chain} k={k}: |C_k|={n}, |C_(k-1)|={m}; #eig>=9/8: {above}; #eig in [9/8,23/20]: {inwin}; "
                  f"max bulk {bulk:.4f}, min band {top:.4f}; lmin(Q1)={eq1:.4f} (>= {float(RHO_L):.4f}); "
                  f"lmin(20M)={em:.2f} (>= {float(RHO_U):.2f}); global residual margins {float(ml):.4f}, {float(mu):.3f}  "
                  f"{'OK' if ok else 'FAIL'}  [{time.time() - t0:.0f}s]", flush=True)
            assert ok
    print("cross-checks passed")


if __name__ == '__main__':
    main()
```

### A.12 `tests/census.py`, `tests/sizes.py`, `tests/equivariance.py`, `tests/boundary_seeds.py`, `tests/attachments.py`

```python
"""Census of the finite 'interior' clusters of the level-K patch: components of the same-species nearest-neighbour graph
all of whose sites have all six lattice neighbours present in the patch.  Such a component is a genuine cluster of the
infinite tiling (no later site can occupy a position that is already occupied).  The sizes found are exactly those of the
two chains.   usage: python -m tests.census [K]"""
import collections
import sys

from ww.core import D, build_patch

K = int(sys.argv[1]) if len(sys.argv) > 1 else 9
patch = build_patch(K)
seen = set()
sizes = collections.Counter()
for p in patch:
    if p in seen:
        continue
    t = patch[p][1]
    comp, stack = [p], [p]
    seen.add(p)
    while stack:
        x = stack.pop()
        for d in D:
            y = (x[0] + d[0], x[1] + d[1])
            if y in patch and y not in seen and patch[y][1] == t:
                seen.add(y)
                stack.append(y)
                comp.append(y)
    if all((s[0] + d[0], s[1] + d[1]) in patch for s in comp for d in D):
        sizes[len(comp)] += 1
A = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605]
B = [2, 10, 42, 165, 632, 2396, 9049, 34135]
print(f"patch level {K}: {len(patch)} sites; interior finite clusters by size (size: number):")
print(dict(sorted(sizes.items())))
extra = set(sizes) - set(A) - set(B)
print("sizes outside the two chains A, B:", sorted(extra))
assert not extra
```

```python
"""Brute-force sizes |C_k| of both chains (numpy inflation of the actual clusters), for comparison with the exact
substitution of ww.counts.   usage: python -m tests.sizes [LA LB]"""
import sys
from tests.cluster_numpy import chain_levels

LA = int(sys.argv[1]) if len(sys.argv) > 1 else 9
LB = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for chain, L in (('A', LA), ('B', LB)):
    sizes = []
    for k, prev, cur in chain_levels(chain, L):
        if k == 1:
            sizes.append(int(prev[4].sum()))
        sizes.append(int(cur[4].sum()))
    print(f"chain {chain}: |C_0..C_{L}| =", sizes, flush=True)
```

```python
"""The children map on windows commutes with the group G = Z_6 x Z_2 (needed to work modulo G)."""
import random
import numpy as np
from ww.core import Frame, geometric_check
from ww.closure import all_sites_closure, child_windows_plain

F = Frame(7)
assert geometric_check(7)
plain = all_sites_closure(F)
random.seed(0)
bad = 0
for v in random.sample(plain, 200):
    ch = [F.from_dict(w) for w in child_windows_plain(F, F.to_dict(v))]
    vr = v
    for k in range(1, 6):
        vr = F.rotate(vr)
        for s in (0, 1):
            vv = F.flip(vr) if s else vr
            for a, w in zip(ch, child_windows_plain(F, F.to_dict(vv))):
                img = a
                for _ in range(k):
                    img = F.rotate(img)
                img = F.flip(img) if s else img
                bad += not np.array_equal(img, F.from_dict(w))
print("windows tested: 200 x 11 group elements x 4 children; equivariance violations:", bad)
assert bad == 0
```

```python
"""Chains that start next to a hole of the tiling (singleton clusters with a lattice neighbour that is not a site): sizes
|C_k| from the substitution of ww.counts, using the same cells as ww.verify --all-seeds.  These chains are "on the boundary"
and are not part of the statement; they are listed only to show how the interior chains differ from them.
usage: python -m tests.boundary_seeds"""
import collections

from ww.closure import cells as build_cells
from ww.counts import children_types, sizes

F, plain, seeds, wins = build_cells(7, interior_only=False)
_, _, seeds_in, _ = build_cells(7, interior_only=True)
interior = {v.tobytes() for _, v in seeds_in}
type_index = {v.tobytes(): i for i, v in enumerate(wins)}
trans = [children_types(F, F.to_dict(v), type_index) for v in wins]
seq = collections.Counter()
for c, v in seeds:
    if c == 'A' and v.tobytes() not in interior:
        seq[tuple(sizes(trans, {type_index[v.tobytes()]: 1}, 7))] += 1
print("singleton seed windows adjacent to a hole:", sum(seq.values()))
for s, n in sorted(seq.items(), reverse=True):
    print(f"  {n:3d} window classes: |C_0..C_7| = {list(s)}")
```

```python
"""How the extras (flipped fourth children) are attached to the rest of the actual clusters.  The statement of the problem
says "attached as a matching"; the proof of FULL_PROOF.md does not use that phrase (only H0-H2).  For every extra of C_k the
number of lattice neighbours of the same species that are triangle vertices, and the number that are extras.
usage: python -m tests.attachments"""
import collections

from tests.cluster_numpy import chain_levels
from tests.crosscheck import cluster_data

for chain, L in (('A', 4), ('B', 3)):
    for k, prev, cur in chain_levels(chain, L):
        cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
        stat = collections.Counter()
        for i, p in enumerate(cs):
            if child[p] == 3:
                nt = sum(1 for j in adj[i] if child[cs[j]] != 3)
                ne = sum(1 for j in adj[i] if child[cs[j]] == 3)
                stat[(nt, ne)] += 1
        print(f"chain {chain} k={k}: |C_k| = {len(cs)}; extras by (#triangle-vertex neighbours, #extra neighbours): "
              f"{dict(sorted(stat.items()))}")
```

---

## Appendix B — program outputs (verbatim)

### B.1 `python -m ww.verify` (the proof: C1–C3 for the two chains with all neighbours present)

```text
all windows of the tiling (radius 7, mod G): 3542
seeds: 103 singleton windows (chain A), 22 pair windows (chain B)
inflation cells (closure under inflation, mod G): 3223   [8s]
lower edge: 896 pieces, all symmetric and positive definite (exact), smallest pivot 1e-06
upper edge: 1028 pieces, all symmetric and positive definite (exact), smallest pivot 0.0102
  cell 0/3223  [10s]
  cell 500/3223  [13s]
  cell 1000/3223  [17s]
  cell 1500/3223  [21s]
  cell 2000/3223  [26s]
  cell 2500/3223  [30s]
  cell 3000/3223  [33s]
lower edge: 18473 rows checked; min row margin rho_L = 85844911/1073741824 = 0.079949 > 0
upper edge: 3223 rows checked; min row margin rho_U = 696063391/8388608 = 82.977222 > 0 (units of 20*M)
=> Q1 = 9/8 - A + 2 sum_v (A u_v)(A u_v)^T  >=  0.0799 * I
=> M = U^T p5(A) U  >=  4.1489 * I     (u_v = indicator of the triangle of v, not normalised)
VERIFIED in 35s
```

### B.2 `python -m ww.verify --all-seeds` (also the chains that start next to holes)

```text
all windows of the tiling (radius 7, mod G): 3542
seeds: 164 singleton windows (chain A), 22 pair windows (chain B)
inflation cells (closure under inflation, mod G): 3440   [9s]
lower edge: 896 pieces, all symmetric and positive definite (exact), smallest pivot 1e-06
upper edge: 1028 pieces, all symmetric and positive definite (exact), smallest pivot 0.0102
  cell 0/3440  [11s]
  cell 500/3440  [13s]
  cell 1000/3440  [17s]
  cell 1500/3440  [21s]
  cell 2000/3440  [25s]
  cell 2500/3440  [29s]
  cell 3000/3440  [33s]
lower edge: 19691 rows checked; min row margin rho_L = 85844911/1073741824 = 0.079949 > 0
upper edge: 3440 rows checked; min row margin rho_U = 696063391/8388608 = 82.977222 > 0 (units of 20*M)
=> Q1 = 9/8 - A + 2 sum_v (A u_v)(A u_v)^T  >=  0.0799 * I
=> M = U^T p5(A) U  >=  4.1489 * I     (u_v = indicator of the triangle of v, not normalised)
VERIFIED in 37s
```

### B.3 `python -m ww.counts` (sizes, recurrence, Perron root, labels)

```text
exact Krylov dimension of the substitution started from a chain-A seed (all 103 seeds): [12]  (=> |C_k| obeys a linear recurrence of order <= 12 for ALL k; the minimal one has order 6)
exact Krylov dimension of the substitution started from a chain-B seed (all matched pairs): [13]
chain A: |C_k| = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605] ...   (satisfies the recurrence with characteristic polynomial x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x-1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2))
chain B: |C_k| = [2, 10, 42, 165, 632, 2396, 9049, 34135, 128737, 485549] ...   (satisfies the recurrence with characteristic polynomial x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x-1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2))
Perron root theta = 3.772538098857,  lim |C_(k-1)|/|C_k| = 1/theta = 0.265073532406,  1 - 1/theta = 0.734926467594
finite-k labels 1 - |C_(k-1)|/|C_k|  (A): ['3/4', '3/4', '47/63', '26/35', '696/941', '2642/3583', '9991/13574']
finite-k labels 1 - |C_(k-1)|/|C_k|  (B): ['4/5', '16/21', '41/55', '467/632', '441/599', '6653/9049', '25086/34135']
quintic irreducible over Q: True
```

### B.4 `python -m tests.crosscheck 6 5` (brute-force clusters: eigenvalues, forms, global exact residuals)

```text
chain A k=1: |C_k|=4, |C_(k-1)|=1; #eig>=9/8: 1; #eig in [9/8,23/20]: 0; max bulk 0.3111, min band 2.1701; lmin(Q1)=0.8104 (>= 0.0799); lmin(20M)=20411.00 (>= 82.98); global residual margins 0.0803, 83.252  OK  [0s]
chain A k=2: |C_k|=16, |C_(k-1)|=4; #eig>=9/8: 4; #eig in [9/8,23/20]: 0; max bulk 1.0000, min band 1.6385; lmin(Q1)=0.0969 (>= 0.0799); lmin(20M)=4492.70 (>= 82.98); global residual margins 0.0800, 83.222  OK  [0s]
chain A k=3: |C_k|=63, |C_(k-1)|=16; #eig>=9/8: 16; #eig in [9/8,23/20]: 0; max bulk 1.0000, min band 1.6097; lmin(Q1)=0.0874 (>= 0.0799); lmin(20M)=3985.17 (>= 82.98); global residual margins 0.0799, 83.182  OK  [0s]
chain A k=4: |C_k|=245, |C_(k-1)|=63; #eig>=9/8: 63; #eig in [9/8,23/20]: 0; max bulk 1.0000, min band 1.4427; lmin(Q1)=0.0871 (>= 0.0799); lmin(20M)=810.45 (>= 82.98); global residual margins 0.0799, 83.162  OK  [0s]
chain A k=5: |C_k|=941, |C_(k-1)|=245; #eig>=9/8: 245; #eig in [9/8,23/20]: 0; max bulk 1.0070, min band 1.3459; lmin(Q1)=0.0803 (>= 0.0799); lmin(20M)=254.58 (>= 82.98); global residual margins 0.0799, 83.098  OK  [1s]
chain A k=6: |C_k|=3583, |C_(k-1)|=941; #eig>=9/8: 941; #eig in [9/8,23/20]: 0; max bulk 1.0070, min band 1.2884; lmin(Q1)=0.0803 (>= 0.0799); lmin(20M)=125.34 (>= 82.98); global residual margins 0.0799, 83.078  OK  [17s]
chain B k=1: |C_k|=10, |C_(k-1)|=2; #eig>=9/8: 2; #eig in [9/8,23/20]: 0; max bulk 0.8575, min band 2.1656; lmin(Q1)=0.1794 (>= 0.0799); lmin(20M)=13149.80 (>= 82.98); global residual margins 0.0803, 83.242  OK  [0s]
chain B k=2: |C_k|=42, |C_(k-1)|=10; #eig>=9/8: 10; #eig in [9/8,23/20]: 0; max bulk 1.0000, min band 1.5197; lmin(Q1)=0.0964 (>= 0.0799); lmin(20M)=1040.05 (>= 82.98); global residual margins 0.0803, 83.172  OK  [0s]
chain B k=3: |C_k|=165, |C_(k-1)|=42; #eig>=9/8: 42; #eig in [9/8,23/20]: 0; max bulk 1.0070, min band 1.3510; lmin(Q1)=0.0803 (>= 0.0799); lmin(20M)=259.23 (>= 82.98); global residual margins 0.0802, 83.116  OK  [0s]
chain B k=4: |C_k|=632, |C_(k-1)|=165; #eig>=9/8: 165; #eig in [9/8,23/20]: 0; max bulk 1.0070, min band 1.2884; lmin(Q1)=0.0803 (>= 0.0799); lmin(20M)=125.35 (>= 82.98); global residual margins 0.0799, 83.078  OK  [0s]
chain B k=5: |C_k|=2396, |C_(k-1)|=632; #eig>=9/8: 632; #eig in [9/8,23/20]: 0; max bulk 1.0070, min band 1.2549; lmin(Q1)=0.0803 (>= 0.0799); lmin(20M)=102.90 (>= 82.98); global residual margins 0.0799, 82.977  OK  [5s]
cross-checks passed
```

### B.5 `python -m tests.equivariance`

```text
windows tested: 200 x 11 group elements x 4 children; equivariance violations: 0
```

### B.6 `python -m tests.census 9`

```text
patch level 9: 262144 sites; interior finite clusters by size (size: number):
{1: 8641, 2: 466, 4: 2045, 10: 84, 16: 473, 42: 12, 63: 106, 165: 1, 245: 23, 941: 5, 3583: 1}
sizes outside the two chains A, B: []
```

### B.7 `python -m tests.crosscheck 9 8 structure` (hypotheses H1, H2 on the actual clusters)

```text
chain A k=1: |C_k|=4, |C_(k-1)|=1: H1, H2 hold
chain A k=2: |C_k|=16, |C_(k-1)|=4: H1, H2 hold
chain A k=3: |C_k|=63, |C_(k-1)|=16: H1, H2 hold
chain A k=4: |C_k|=245, |C_(k-1)|=63: H1, H2 hold
chain A k=5: |C_k|=941, |C_(k-1)|=245: H1, H2 hold
chain A k=6: |C_k|=3583, |C_(k-1)|=941: H1, H2 hold
chain A k=7: |C_k|=13574, |C_(k-1)|=3583: H1, H2 hold
chain A k=8: |C_k|=51293, |C_(k-1)|=13574: H1, H2 hold
chain A k=9: |C_k|=193605, |C_(k-1)|=51293: H1, H2 hold
chain B k=1: |C_k|=10, |C_(k-1)|=2: H1, H2 hold
chain B k=2: |C_k|=42, |C_(k-1)|=10: H1, H2 hold
chain B k=3: |C_k|=165, |C_(k-1)|=42: H1, H2 hold
chain B k=4: |C_k|=632, |C_(k-1)|=165: H1, H2 hold
chain B k=5: |C_k|=2396, |C_(k-1)|=632: H1, H2 hold
chain B k=6: |C_k|=9049, |C_(k-1)|=2396: H1, H2 hold
chain B k=7: |C_k|=34135, |C_(k-1)|=9049: H1, H2 hold
chain B k=8: |C_k|=128737, |C_(k-1)|=34135: H1, H2 hold
```

### B.8 `python -m tools.cert_stats`

```text
certificate L: type radius Rt = 2, entries are integers over 2^30, 896 pieces
  piece sizes (size: number of pieces): {1: 37, 2: 6, 4: 158, 5: 17, 6: 100, 7: 65, 8: 67, 9: 202, 10: 90, 11: 41, 12: 59, 13: 10, 14: 35, 15: 9}
  pieces by child index j of the centre: {0: 225, 1: 225, 2: 225, 3: 221}
  largest |entry| = 4.51129
  sha256 97fb3ccaf445647249d387f5cfb80dbba8977129b0bdc9c1a5c939d42061a949  (gzip file, 248300 bytes)
  example smallest piece: j = 3, 1x1, offsets [[0, 0]]
      [1.0447] 
  example largest piece: j = 0, 15x15, offsets [[-2, 0], [-2, 2], [-1, 0], [-1, 1]]...
      [0.66737, 0.59829, 1.00955, 0.87467, 0.58502, 0.15757] ...
      [0.59829, 0.88382, 0.93122, 1.33084, 0.96534, 0.07445] ...
      [1.00955, 0.93122, 1.69789, 1.45931, 0.84544, 0.55175] ...
certificate U: type radius Rt = 4, entries are integers over 2^24, 1028 pieces
  piece sizes (size: number of pieces): {1: 64, 2: 8, 3: 50, 4: 169, 5: 121, 6: 73, 7: 151, 8: 134, 9: 167, 10: 25, 11: 30, 12: 28, 13: 8}
  largest |entry| = 89218.3
  sha256 3f9d974bd7b5f9422520de28a31192db7592864e908053bf0f6c8ddf0cf9955c  (gzip file, 242782 bytes)
  example smallest piece: j = -1, 1x1, offsets [[0, 0]]
      [20327.74766] 
  example largest piece: j = -1, 13x13, offsets [[-2, 0], [-1, 0], [-1, 2], [0, 0]]...
      [1041.1683, 4154.78752, 307.56205, 415.98933, -128.3817, 81.60876] ...
      [4154.78752, 18618.52521, 2584.75364, 4389.06441, 2905.14042, 1848.46554] ...
      [307.56205, 2584.75364, 22037.62219, 12900.68658, 29359.92954, 16734.26716] ...
```

### B.9 `python -m tests.sizes 9 8`

```text
chain A: |C_0..C_9| = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605]
chain B: |C_0..C_8| = [2, 10, 42, 165, 632, 2396, 9049, 34135, 128737]
```

### B.10 `python -m tools.algebra`

```text
x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x - 1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2): verified
f mod 3 = [2, 0, 2, 0, 2, 1] (ascending coefficients)
f mod 3 has no monic factor of degree 1 or 2 (all 3 + 9 candidates tested): f is irreducible over F_3, hence over Q
Sturm: f has exactly 3 real roots
f at -1/2, -2/5, 39/100, 2/5, 3, 4: ['-43/32', '1408/3125', '458665499/10000000000', '-28/3125', '-16', '22']
non-real roots: |z|^2 = 2/(theta |r1| r2) < 2/(3 * 2/5 * 39/100) = 500/117 = 4.2735 < 9, so |z| < 3 < theta
theta        = 3.7725380988570505404160270370041892583877
1/theta      = 0.2650735324059326665550111967750136168646
1 - 1/theta  = 0.7349264675940673334449888032249863831354
y = 1/theta is a root of (1)y^0 + (-7)y^1 + (15)y^2 + (-10)y^3 + (-3)y^4 + (2)y^5 = 0
l = 1 - 1/theta is a root of (-2)l^0 + (9)l^1 + (-13)l^2 + (2)l^3 + (7)l^4 + (-2)l^5 = 0
chain A: the recurrence started from s_0..s_5 reproduces the brute-force sizes s_6..s_9
   ratios s_(k-1)/s_k: ['0.25000', '0.25000', '0.25397', '0.25714', '0.26036', '0.26263', '0.26396', '0.26464', '0.26494']
chain B: the recurrence started from s_0..s_5 reproduces the brute-force sizes s_6..s_8
   ratios s_(k-1)/s_k: ['0.20000', '0.23810', '0.25455', '0.26108', '0.26377', '0.26478', '0.26509', '0.26515']
```

### B.11 `python -m tests.boundary_seeds`

```text
singleton seed windows adjacent to a hole: 61
    9 window classes: |C_0..C_7| = [1, 4, 16, 62, 236, 889, 3336, 12515]
   36 window classes: |C_0..C_7| = [1, 4, 15, 56, 210, 789, 2966, 11152]
   16 window classes: |C_0..C_7| = [1, 4, 15, 55, 201, 737, 2719, 10093]
```

### B.12 `python -m tests.attachments`

```text
chain A k=1: |C_k| = 4; extras by (#triangle-vertex neighbours, #extra neighbours): {(1, 0): 1}
chain A k=2: |C_k| = 16; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 1, (1, 0): 2, (2, 1): 1}
chain A k=3: |C_k| = 63; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 3, (1, 0): 4, (2, 0): 1, (2, 1): 3, (3, 0): 4}
chain A k=4: |C_k| = 245; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 10, (1, 0): 12, (2, 0): 7, (2, 1): 10, (3, 0): 17}
chain B k=1: |C_k| = 10; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 1, (1, 0): 1, (2, 1): 1, (3, 0): 1}
chain B k=2: |C_k| = 42; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 3, (1, 0): 3, (2, 0): 1, (2, 1): 3, (3, 0): 2}
chain B k=3: |C_k| = 165; extras by (#triangle-vertex neighbours, #extra neighbours): {(0, 1): 8, (1, 0): 9, (2, 0): 5, (2, 1): 8, (3, 0): 9}
```

---

## Appendix C — fingerprints (SHA-256)

| file | sha256 |
|------|--------|
| `data/cert_L.json.gz` (896 pieces for Q₁) | 97fb3ccaf445647249d387f5cfb80dbba8977129b0bdc9c1a5c939d42061a949 |
| `data/cert_U.json.gz` (1028 pieces for 20M) | 3f9d974bd7b5f9422520de28a31192db7592864e908053bf0f6c8ddf0cf9955c |
| `ww/core.py` | b4b70dd157a7a138d1719e56281992aa1605d35f74ed18781a56b2e502e6f4c1 |
| `ww/closure.py` | 8ee6c028508d36e592bed4634c22e0640cc933366e7a8a2a5848858ce681ac5f |
| `ww/cells.py` | 9837fd81cd261dfc4419a00fa641818ddf7f65106473d2b29972a467d955e9f8 |
| `ww/verify.py` | 1b872b78cc9264e36131bafdf1a5f4c1f7bc9f92c4cd2b89bd214f90053318de |
| `ww/counts.py` | f98e8d7e205470725890273a22d9c70afa17d0da46d8a9259439a0a1633fee01 |
| `tools/algebra.py` | aac165e66330a24e76da66d830302ddee5a923c9bd017e67e8c28fab0d4dfc1e |
| `tools/gen_constraints.py` | bfde3c83396a4e8341002b6b5b2e3afb161a19a4ab4c20eda555f8b727cadcc6 |
| `tools/solve_sdp.py` | a398b5a80f49426cf457d6f362c0cbc60f978a23d53b037e09f8e9261bb8fbdd |
| `tools/cert_stats.py` | 4d0a6e23f17471e640e9bf974a3e30ec92f6657371f9ea78fa42be78f37d38f0 |
| `tests/cluster_numpy.py` | f038c56c28436f80d4bd826c654fdba113865e17425a3c13a7d075c9bd5149c6 |
| `tests/crosscheck.py` | a731ff1e4dd5d91decf4762bf3b226f72afac6ee143b4ee4226b525c14cf4adc |
| `tests/census.py` | 20de0b102bfffdf48f5f07cae4d26233e937ed14eaa83f0b58a4d42483f3ced0 |
| `tests/sizes.py` | ffb67f980c3e2113b5d4aab3d6a51d4aa6102295ed4b5c4bbabcebf671668e97 |
| `tests/equivariance.py` | ccab8d28506df8e4fcd7fb7cf06dc5c31d60ce2c1b706aa2191ab922f2f26d1c |
| `tests/boundary_seeds.py` | f6a19e7a86aa04cb60ac8ee94df7336597cccf473da8d59fdcaff21f18f832ac |
| `tests/attachments.py` | 356f313e0ceb69cd737812e4537e84e9fdcc24c86733671dc96b0846768c25da |

**The certificate data are the only part of the computation that is not inlined in this file.** The two files are gzip-compressed JSON of about 0.25 MB each (about 63 000 resp. 50 000 stored integers, 35 000 resp. 28 000 distinct by symmetry): {"form", "Rt", "bits", "types": [{"j": child index (−1 for the upper edge), "key": hex of the canonical radius-Rt window, "offs": canonical offsets, "N": symmetric integer matrix whose entries are multiples of 2^−bits}]}. A reader who wants to re-run the verification needs them (repository directory `data/`, the hashes above identify them); alternatively the pieces can be regenerated from scratch with Appendix A.7–A.8 (the numerical solution differs slightly from run to run — what matters is only that the exact check of Appendix A.4 passes). Statistics and example pieces are in B.8.

## Appendix D — how to reproduce

Software used for the outputs in Appendix B: Python 3.13.17, numpy 2.5.3, sympy 1.14.0, scipy 1.18.1, cvxpy 1.9.3, scs 3.3.1 (the last three only for regenerating certificates).

```
python -m ww.verify                  # C1–C3, ≈ 35 s                    (output B.1)
python -m ww.verify --all-seeds      # also the chains next to holes    (output B.2)
python -m ww.counts                  # §7, ≈ 1 min                      (output B.3)
python -m tests.crosscheck 6 5       # §6                               (output B.4)
python -m tests.equivariance         #                                  (output B.5)
python -m tests.census 9             #                                  (output B.6)
python -m tests.crosscheck 9 8 structure                                (output B.7)
python -m tools.cert_stats           #                                  (output B.8)
python -m tests.sizes 9 8            #                                  (output B.9)
python -m tools.algebra              #                                  (output B.10)
python -m tests.boundary_seeds       #                                  (output B.11)
python -m tests.attachments          #                                  (output B.12)

# regenerate the certificates (optional; needs scipy, cvxpy, scs):
python -m tools.gen_constraints L 2 L.pkl --all-seeds
python -m tools.gen_constraints U 4 U.pkl --all-seeds
python -m tools.solve_sdp L.pkl cert_L.json.gz L --tol 1e-5 --bits 30 --margin 1e-6
python -m tools.solve_sdp U.pkl cert_U.json.gz U --tol 1e-6 --scale 1000 --bits 24 --margin 1e-2
python -m ww.verify --all-seeds --data <directory holding the two new files>
```
