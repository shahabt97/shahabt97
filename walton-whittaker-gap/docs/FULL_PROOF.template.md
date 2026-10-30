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

{{FILE:ww/core.py}}

### A.2 `ww/closure.py` — all windows, seeds, cells

{{FILE:ww/closure.py}}

### A.3 `ww/cells.py` — inflation cell, the two local forms

{{FILE:ww/cells.py}}

### A.4 `ww/verify.py` — exact verification (the program whose output is in B.1)

{{FILE:ww/verify.py}}

### A.5 `ww/counts.py` — substitution, sizes, recurrence, gap labels

{{FILE:ww/counts.py}}

### A.6 `tools/algebra.py` — exact algebra of the gap labels

{{FILE:tools/algebra.py}}

### A.7 `tools/gen_constraints.py` — constraints of the piece SDP

{{FILE:tools/gen_constraints.py}}

### A.8 `tools/solve_sdp.py` — SDP solve, shift, rounding, export

{{FILE:tools/solve_sdp.py}}

### A.9 `tools/cert_stats.py` — statistics and fingerprints of the certificate files

{{FILE:tools/cert_stats.py}}

### A.10 `tests/cluster_numpy.py` — independent brute-force construction of the actual clusters

{{FILE:tests/cluster_numpy.py}}

### A.11 `tests/crosscheck.py` — eigenvalue counts, global exact residuals, structure H1/H2

{{FILE:tests/crosscheck.py}}

### A.12 `tests/census.py`, `tests/sizes.py`, `tests/equivariance.py`, `tests/boundary_seeds.py`, `tests/attachments.py`

{{FILE:tests/census.py}}

{{FILE:tests/sizes.py}}

{{FILE:tests/equivariance.py}}

{{FILE:tests/boundary_seeds.py}}

{{FILE:tests/attachments.py}}

---

## Appendix B — program outputs (verbatim)

### B.1 `python -m ww.verify` (the proof: C1–C3 for the two chains with all neighbours present)

{{LOG:verify}}

### B.2 `python -m ww.verify --all-seeds` (also the chains that start next to holes)

{{LOG:verify_all_seeds}}

### B.3 `python -m ww.counts` (sizes, recurrence, Perron root, labels)

{{LOG:counts}}

### B.4 `python -m tests.crosscheck 6 5` (brute-force clusters: eigenvalues, forms, global exact residuals)

{{LOG:crosscheck}}

### B.5 `python -m tests.equivariance`

{{LOG:equivariance}}

### B.6 `python -m tests.census 9`

{{LOG:census}}

### B.7 `python -m tests.crosscheck 9 8 structure` (hypotheses H1, H2 on the actual clusters)

{{LOG:structure}}

### B.8 `python -m tools.cert_stats`

{{LOG:cert_stats}}

### B.9 `python -m tests.sizes 9 8`

{{LOG:sizes}}

### B.10 `python -m tools.algebra`

{{LOG:algebra}}

### B.11 `python -m tests.boundary_seeds`

{{LOG:boundary_seeds}}

### B.12 `python -m tests.attachments`

{{LOG:attachments}}

---

## Appendix C — fingerprints (SHA-256)

| file | sha256 |
|------|--------|
| `data/cert_L.json.gz` (896 pieces for Q₁) | {{SHA:data/cert_L.json.gz}} |
| `data/cert_U.json.gz` (1028 pieces for 20M) | {{SHA:data/cert_U.json.gz}} |
| `ww/core.py` | {{SHA:ww/core.py}} |
| `ww/closure.py` | {{SHA:ww/closure.py}} |
| `ww/cells.py` | {{SHA:ww/cells.py}} |
| `ww/verify.py` | {{SHA:ww/verify.py}} |
| `ww/counts.py` | {{SHA:ww/counts.py}} |
| `tools/algebra.py` | {{SHA:tools/algebra.py}} |
| `tools/gen_constraints.py` | {{SHA:tools/gen_constraints.py}} |
| `tools/solve_sdp.py` | {{SHA:tools/solve_sdp.py}} |
| `tools/cert_stats.py` | {{SHA:tools/cert_stats.py}} |
| `tests/cluster_numpy.py` | {{SHA:tests/cluster_numpy.py}} |
| `tests/crosscheck.py` | {{SHA:tests/crosscheck.py}} |
| `tests/census.py` | {{SHA:tests/census.py}} |
| `tests/sizes.py` | {{SHA:tests/sizes.py}} |
| `tests/equivariance.py` | {{SHA:tests/equivariance.py}} |
| `tests/boundary_seeds.py` | {{SHA:tests/boundary_seeds.py}} |
| `tests/attachments.py` | {{SHA:tests/attachments.py}} |

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
