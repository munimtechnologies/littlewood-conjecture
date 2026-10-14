# The Littlewood conjecture in October 2026

**What is proved, what changed in 2026, and what is still open about simultaneous approximation
by a product of distances.** This is an exposition with complete proofs of the elementary reductions
(the facts that leave only badly approximable, linearly independent pairs) and of the step that turns
a counterexample into a bounded orbit, plus a small, tested Python companion (`littlewood.py`) that
computes certified continued fractions and certified Littlewood products.

By MV2 of Munim, Inc. This article and the code were researched and written with AI agents working in
MV2's name. Every statement below is either proved here or cited to a source. Section 10 separates
what is proved from what is open, and section 11 lists what is new here. Corrections are welcome by
issue or pull request.

Version 1.0, 2026-10-14. Text CC BY 4.0, code Apache-2.0.

Related problem records: Justin Sun Prize JSP-000016 (MathDB #316099, "Littlewood conjecture"),
Wikipedia "Littlewood conjecture".

---

## 1. The question

Write ‖t‖ for the distance from a real number t to the nearest integer. **Littlewood's conjecture**
(posed by J. E. Littlewood around 1930) says:

> For all real numbers α and β,  liminf_{q→∞} q · ‖qα‖ · ‖qβ‖ = 0.

Dirichlet's theorem gives, for one number, infinitely many q with q‖qα‖ < 1, and for two numbers
infinitely many q with ‖qα‖, ‖qβ‖ < q^{−1/2}. Both give a product q‖qα‖‖qβ‖ that is merely bounded.
The conjecture asks for more: the **product** of the two errors, scaled by q, can be made as small as
we like. Two equivalent forms are used below.

- **The infimum form.** Put m(α, β) = inf_{q ≥ 1} q‖qα‖‖qβ‖. The conjecture is equivalent to
  m(α, β) = 0 for all pairs (Lemma 1), i.e. to m_LC = sup_{α,β} m(α, β) being 0.
- **Products of linear forms.** m(α, β) = inf |x (xα − y)(xβ − z)| over integers x ≠ 0, y, z (take y,
  z nearest to xα, xβ). So the conjecture asks whether this product of three linear forms in
  (x, y, z) takes values arbitrarily close to 0 at integer points with x ≠ 0.

## 2. Elementary reductions, with proofs

**Lemma 1 (infimum versus liminf).** liminf q‖qα‖‖qβ‖ = 0 if and only if m(α, β) = 0.
*Proof.* One direction is clear. Conversely, if some term is 0 then α or β is rational, and then
infinitely many terms are 0 (multiples of the denominator). If no term is 0, all terms are positive,
so an infimum of 0 is not attained and is approached by infinitely many q. ∎

**Continued fractions.** Let α = [a₀; a₁, a₂, ...] be irrational, with convergents pₖ/qₖ
(q₋₁ = 0, q₀ = 1, qₖ₊₁ = aₖ₊₁qₖ + qₖ₋₁). Writing α = [a₀; ..., aₖ, xₖ₊₁] with complete quotient
xₖ₊₁ ∈ (aₖ₊₁, aₖ₊₁ + 1), the identity qₖα − pₖ = ±1/(xₖ₊₁qₖ + qₖ₋₁) gives, for k ≥ 1,

    1 / (qₖ₊₁ + qₖ)  <  ‖qₖα‖  <  1 / qₖ₊₁.                    (Legendre's bounds)

α is **badly approximable** (α ∈ **Bad**) if inf_q q‖qα‖ > 0. By the bounds,
qₖ‖qₖα‖ < qₖ/qₖ₊₁ < 1/aₖ₊₁, so α ∉ Bad when the partial quotients are unbounded. (Conversely, bounded
partial quotients imply α ∈ Bad, because the best approximations are convergents; this is classical.)

**Lemma 2 (one coordinate suffices).** If α is rational or α ∉ Bad, the conjecture holds for (α, β)
for every β.
*Proof.* Rational α: shown in Lemma 1. Otherwise inf_q q‖qα‖ = 0 and ‖qβ‖ ≤ 1/2, so
m(α, β) ≤ ½ inf_q q‖qα‖ = 0; apply Lemma 1. ∎

**Lemma 3 (linear dependence, with an explicit rate).** If a + bα + cβ = 0 for integers a, b, c not all
zero, the conjecture holds for (α, β). More precisely, if c ≠ 0 and α is irrational, then for every
convergent denominator qₖ of α (k ≥ 1) the integer n = |c|qₖ satisfies

    n ‖nα‖ ‖nβ‖  ≤  c² |b| qₖ / qₖ₊₁²  <  c² |b| / qₖ₊₁.

*Proof.* If c = 0 then α is rational (b ≠ 0), and Lemma 2 applies; if α is rational, likewise.
Otherwise ‖nα‖ = ‖|c|·qₖα‖ ≤ |c|·‖qₖα‖ < |c|/qₖ₊₁, using ‖mt‖ ≤ |m|·‖t‖ for integers m. And
nβ = ±qₖ(cβ) = ∓qₖ(a + bα), so ‖nβ‖ = ‖qₖbα‖ ≤ |b|·‖qₖα‖ < |b|/qₖ₊₁. Multiply by n = |c|qₖ. ∎

The test suite checks this for β = (1 + 2√2)/3, α = √2: with qₖ = 2378, qₖ₊₁ = 5741 and n = 7134 the
bound is 18 · 2378 / 5741² < 0.0013, and the computed value is about 0.00095.

**Lemma 4 (a universal bound).** m(α, β) < 1/(2√5) = 0.2236... for every pair.
*Proof.* By Hurwitz's theorem every irrational α has infinitely many q with q‖qα‖ < 1/√5 (Borel: one of
any three consecutive convergents works); multiply by ‖qβ‖ ≤ 1/2. Rational α gives m = 0. ∎

The real question is how small m(α, β) can be made *uniformly*: Badziahin (2016) proved m_LC < 1/19 and
Matala-Aho, Törmä and Vapa (2026) proved m_LC < 1/26 (section 7). A proof that m_LC = 0 is the conjecture.

**Lemma 5 (quadratic irrationals are badly approximable).** If α is a real root of an irreducible
f(x) = Ax² + Bx + C ∈ ℤ[x], then α ∈ Bad.
*Proof.* For integers p, q ≥ 1, q²f(p/q) is a nonzero integer, so |f(p/q)| ≥ 1/q². If |α − p/q| ≤ 1,
then |f(p/q)| = |A|·|p/q − α|·|p/q − α'| ≤ |A|(|α − α'| + 1)·|α − p/q|. Hence
q²|α − p/q| ≥ 1/(|A|(|α − α'| + 1)) for all p/q close to α, and q‖qα‖ is bounded below. ∎

**Where this leaves the problem.** A counterexample must be a pair of badly approximable numbers α, β
with 1, α, β linearly independent over ℚ. **(√2, √3)** is such a pair (Lemma 5; 1, √2, √3 are
independent), and the conjecture is **not known** for it (Gowers 2009; we found no later result
settling it). By Jarník, Bad has Hausdorff dimension 1, so these reductions alone leave a set of
dimension 2.

## 3. Almost every pair, and the algebraic pairs

- **Almost every pair.** Bad has Lebesgue measure 0 (Khintchine 1924; Borel–Bernstein), so by Lemma 2
  the possible counterexamples form a null set. Gallagher (1962) proved much more: for almost every
  (α, β), liminf q (log q)² ‖qα‖‖qβ‖ = 0. Chow and Yang (Invent. Math. 2024) proved Gallagher-type
  statements on planar lines, using effective equidistribution.
- **Same cubic field.** Cassels and Swinnerton-Dyer (1955) proved the conjecture when 1, α, β span a
  real cubic field, by proving the dual statement inf_{m,n≠0} |mn|·‖mα + nβ‖ = 0 and transferring it.
  Peck (1961) sharpened this: for a basis 1, α₁, ..., α_d of a real number field of degree d + 1,
  liminf q log q ‖qα₁‖...‖qα_d‖ < ∞. In September 2026 Dhanda, Flynn and Haynes gave a geometric proof
  of Peck's theorem from one application of Minkowski's convex body theorem (arXiv 2609.29469), and
  Bugeaud and de Mathan proved refinements toward Peck's conjecture on splitting the logarithmic saving
  unequally between coordinates (arXiv 2609.29360).
- **Two different quadratic fields**, such as (√2, √3), are not covered: they live in the quartic field
  ℚ(√2, √3), and Peck-type results for its natural basis 1, √2, √3, √6 involve the third number √6.
- **A related open problem.** No algebraic number of degree ≥ 3 is known to be badly approximable, and
  none is known not to be. The companion code certifies the first 41 partial quotients of 2^{1/3},
  [1; 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, 2, 1, 3, 4, 1, 1, 2, 14, 3, 12, 1, 15,
  3, 1, 4, 534, 1, 1, 5, 1, 1] (OEIS A002945); the 534 is a famous surprise. If cubic irrationals are
  not badly approximable, as is widely expected, the Cassels–Swinnerton-Dyer theorem follows
  immediately from Lemma 2.

## 4. Dynamics: why the exceptional set is small

Attach to (α, β) the lattice

    Λ = { (x, xα − y, xβ − z) : x, y, z ∈ ℤ }  ⊂ ℝ³

(covolume 1, so a point of X₃ = SL(3, ℝ)/SL(3, ℤ)), and the diagonal matrices
a_{s,t} = diag(e^{−s−t}, e^{s}, e^{t}) with s, t ≥ 0. They preserve the product of the three coordinates.

**Lemma 6 (a counterexample gives a bounded orbit).** If m(α, β) = δ > 0, then every nonzero vector of
every lattice a_{s,t}Λ (s, t ≥ 0) has sup-norm at least min(1, δ^{1/3}). By Mahler's compactness
criterion the orbit {a_{s,t}Λ : s, t ≥ 0} is then bounded in X₃.
*Proof.* Take a nonzero v = a_{s,t}(x, xα − y, xβ − z). If x ≠ 0, the product of the coordinates of v is
x(xα − y)(xβ − z), of absolute value ≥ δ, so its largest coordinate is ≥ δ^{1/3}. If x = 0, then
v = (0, −e^{s}y, −e^{t}z) with (y, z) ≠ (0, 0) and s, t ≥ 0, so ‖v‖ ≥ 1. ∎

So the conjecture would follow from: *no such orbit is bounded*. This is the case d = 3 of
conjectures of Cassels–Swinnerton-Dyer (products of d ≥ 3 linear forms) and Margulis (bounded orbits
of the full diagonal group are compact). Maucourant (2010) showed that the analogous orbit-closure
statement fails for diagonal subgroups of SL(d, ℝ) with d ≥ 6; this does not affect Littlewood.

**Einsiedler, Katok and Lindenstrauss (Annals 2006)** classified the invariant measures of positive
entropy for the diagonal group on X₃ and deduced: **the set of pairs (α, β) for which the conjecture
fails has Hausdorff dimension 0**; in fact it is a countable union of compact sets of box dimension 0.
This is still the strongest general result (Venkatesh's 2008 survey explains the proof). It does not
show the set is empty. Einsiedler and Kleinbock (2007) proved the same dimension statement for the
p-adic version (section 8).

## 5. Explicit pairs and limits of strengthening

- **Pollington–Velani (Acta Math. 2000):** for every α ∈ Bad, there is a set of β ∈ Bad of Hausdorff
  dimension 1 with q log q ‖qα‖‖qβ‖ ≤ 1 for infinitely many q.
- **de Mathan (2003)** and **Adamczewski–Bugeaud (2006):** explicit pairs of badly approximable
  numbers, with 1, α, β independent, that satisfy the conjecture; for each α with bounded partial
  quotients, Adamczewski and Bugeaud construct uncountably many such β.
- **How much can be gained over Littlewood?** Gallagher's (log q)² holds almost everywhere, but
  Bugeaud and Moshchevitin (2011) showed that the pairs with liminf q (log q)² ‖qα‖‖qβ‖ > 0 have full
  Hausdorff dimension, and Badziahin (2013) did the same with q log q log log q. Fregoli and Kleinbock
  (J. Number Theory, 2026) extended Badziahin's result to ℝ^d with (log q)^{d−1} log log q.
- **Inhomogeneous versions.** Shapira (Annals 2011): for almost every (α, β), for all shifts (γ, δ),
  liminf |n|‖nα − γ‖‖nβ − δ‖ = 0; this also holds for pairs from totally real cubic fields. Further
  results on fibres by Chow–Zafeiropoulos (2021), Chow–Technau (2024) and Stefanescu (2025).

## 6. Products of three linear forms in a picture

For a fixed pair, the companion records every q at which q‖qα‖‖qβ‖ reaches a new minimum
(`python littlewood.py records sqrt2 sqrt3 1000000`). A conjecture-sized question becomes concrete:
does this record sequence tend to 0 for (√2, √3)? A computation can only ever show values getting
small; it cannot prove the liminf is 0, and it cannot exhibit a counterexample. One value that can be
checked by hand: at q = 41,

    41√2 = 57.98275...,  41√3 = 71.01408...,   41 · 0.0172439 · 0.0140831 = 0.0099568...

already below the universal bound 1/26 of section 7. (`python littlewood.py value sqrt2 sqrt3 41`
prints a certified interval of width below 10⁻³⁰.)

## 7. What has been computed

- **Universal bounds on m_LC.** Badziahin (Exp. Math. 2016) gave an algorithm that certifies
  inequalities m_LC < ε, and proved m_LC < 1/19 (and an analogous 2-adic bound of 1/9). Matala-Aho,
  Törmä and Vapa (Exp. Math., online 2026-04-02) gave a new algorithm based on the continued fractions
  of α and β and proved **inf_q q‖qα‖‖qβ‖ < 1/26 for all real α, β**. These are upper bounds for the
  worst pair; the conjecture says the worst pair still gives 0.
- **p-adic numerics.** Badziahin's "Estimating lower limit in the p-adic Littlewood conjecture"
  (arXiv 2506.16860, 2025) and Vitorino–Vucusic (Mathematika 2026) give further bounds for the 2-adic
  case.

## 8. The 2025–2026 developments

**(a) The uniform Littlewood conjecture is false.** Bandi, Fregoli and Kleinbock (arXiv 2504.02258,
2025) proposed a uniform strengthening: for every pair,

    lim_{Q→∞}  Q · min_{1≤q≤Q} ‖qα‖‖qβ‖ = 0,

which they knew to hold for almost every pair. Since q ≤ Q, this implies Littlewood's conjecture; it
asks that *every* large Q has a good q ≤ Q, not just infinitely many. In 2026:

- Schleischitz (arXiv 2603.12611, March) **disproved** it: the counterexamples form a residual set (a
  dense G_δ). The disproof is semi-constructive, using effective results toward Zaremba's conjecture
  (Bourgain–Kontorovich). He also disproved a uniform p-adic version.
- Moshchevitin (arXiv 2605.26188, May) gave a short proof of a counterexample using properties of
  Fibonacci numbers.
- Shulga (arXiv 2608.25059, August) proved the counterexamples with badly approximable first coordinate
  form a set of Hausdorff dimension at least 3/2.
- Neckrasov, Wu and Yang (arXiv 2608.24401, August–September) proved the set of counterexamples is
  hyperplane absolute winning, so it has **full** Hausdorff dimension 2.

This does **not** affect the classical conjecture, whose exceptional set has dimension 0 (section 4).
The contrast is the point: the uniform statement fails on a big set; the classical one can fail only on
a set of dimension 0. `python littlewood.py uniform X Y Q` prints the uniform quantity.

**(b) Function fields: the analogy breaks in odd characteristic.** De Mathan and Teulié also proposed
an analogue of the p-adic conjecture (part (c)) for Laurent series over a finite field 𝔽_q, with an
irreducible polynomial P(t) in place of p: the t-adic and P(t)-adic Littlewood conjectures.
Adiceam, Nesharim and Lunnon (Duke 2021) found an explicit counterexample over 𝔽₃. Garrett and
Robertson (Math. Comp. 2026) found computer-assisted counterexamples to the P(t)-adic version over
small fields (characteristic 5, and the characteristic-3 construction in characteristics 7 and 11); Robertson
(Mathematika 2026) gave a transference from t to any irreducible P(t) using number walls. **Lai and
Sprang (arXiv 2606.00633, May 2026) proved that the P(t)-adic Littlewood conjecture fails for every
irreducible P(t) over every ground field of odd characteristic.** Characteristic 2 is open.
Badziahin, Pavlenkov and Zorin (arXiv 2608.22078, August 2026) showed that in the t-adic case a
non-empty exceptional set has infinite 1/log(1/r)-Hausdorff measure, refining the Lai–Sprang
construction for odd q. None of this implies anything about real numbers.

**(c) The p-adic Littlewood conjecture** (de Mathan–Teulié 2004): for every real α and prime p,
inf_q q·|q|_p·‖qα‖ = 0. Proved for quadratic irrationals (de Mathan–Teulié); exceptional set of
dimension 0 (Einsiedler–Kleinbock 2007); still open. Badziahin, Pavlenkov and Zorin (2026) showed that
if the exceptional set is non-empty, it has positive Hausdorff measure for the gauge 1/log(1/r), so it
has the cardinality of the continuum. `python littlewood.py mixed X P Q` records q|q|_P‖qX‖.

## 9. Status as of 2026-10-14

- **Open.** We found no 2025–2026 result proving or disproving the classical conjecture. The most recent
  papers on it (Shulga, August 2026; Neckrasov–Wu–Yang, September 2026) describe it as open, with the
  dimension-0 theorem of Einsiedler–Katok–Lindenstrauss as the best general result; Wikipedia's article
  (edited 2026-10-05) agrees. The Justin Sun Prize MathDB pilot register lists JSP-000016 with no
  solution, Lean or exposition entry. We found no claimed proof on arXiv in 2025–2026.
- **What moved in 2026:** the uniform version was disproved (with full-dimensional counterexamples), the
  P(t)-adic analogue was disproved in every odd characteristic, the universal bound improved from 1/19
  to 1/26, Badziahin's logarithmic obstruction was extended to higher dimensions, and Peck's theorem got
  a geometric proof and refinements.

## 10. Proved versus open, as of 2026-10-14

**Proved:**

1. The conjecture holds when α or β is rational or not badly approximable, and when 1, α, β are linearly
   dependent over ℚ (Lemmas 2 and 3).
2. m(α, β) < 1/26 for every pair (Matala-Aho–Törmä–Vapa 2026; earlier 1/19, Badziahin 2016; 1/(2√5) by
   Hurwitz, Lemma 4).
3. It holds for almost every pair, with the factor (log q)² (Gallagher 1962), and for pairs spanning a
   real cubic field with 1, with the factor log q (Cassels–Swinnerton-Dyer 1955; Peck 1961).
4. The exceptional set has Hausdorff dimension 0 (Einsiedler–Katok–Lindenstrauss 2006); same for the
   p-adic version (Einsiedler–Kleinbock 2007).
5. For each α ∈ Bad, a dimension-1 set of β ∈ Bad satisfies it (Pollington–Velani 2000); explicit
   independent badly approximable pairs satisfy it (de Mathan 2003; Adamczewski–Bugeaud 2006).
6. The uniform Littlewood conjecture is false, with counterexamples of full dimension (2026).
7. The P(t)-adic analogue is false over every field of odd characteristic (Lai–Sprang 2026, building on
   Adiceam–Nesharim–Lunnon 2021 and Garrett–Robertson 2026).

**Open:**

1. Littlewood's conjecture itself: is m(α, β) = 0 for every pair? Equivalently, is m_LC = 0?
2. Any specific independent pair of algebraic numbers from different quadratic fields, starting with
   (√2, √3).
3. The p-adic Littlewood conjecture (de Mathan–Teulié), and the t-adic one in characteristic 2.
4. The Cassels–Swinnerton-Dyer and Margulis conjectures on diagonal orbits in X₃.
5. Whether any algebraic number of degree ≥ 3 is badly approximable.

## 11. What is original here

- Self-contained proofs of the reductions (Lemmas 1–5), including an explicit rate c²|b|qₖ/qₖ₊₁² in the
  linearly dependent case, and of the step from a counterexample to a bounded diagonal orbit (Lemma 6).
- A **certified continued fraction** method from an exact floor oracle: if N·x is available as
  ⌊N x⌋ in integers, the reals whose floor algorithm starts with a₀, ..., a_{j−1} form an interval, so the
  terms shared by the expansions of ⌊2ᵇx⌋/2ᵇ and (⌊2ᵇx⌋ + 1)/2ᵇ are proved terms of x. With integer k-th
  roots this certifies, for example, 41 partial quotients of 2^{1/3}, including the 534.
- **Certified Littlewood products**: since ‖·‖ is 1-Lipschitz and qx lies in [qA/S, (qA + q)/S) for
  A = ⌊Sx⌋, each value q‖qα‖‖qβ‖ comes with a proved interval; with S = 2²⁵⁶ the width is far below
  10⁻³⁰ for q in the millions.
- `littlewood.py` and `test_littlewood.py`: exact reals ((a + b√d)/c and (p/q)^{1/k}), certified
  continued fractions and convergents, certified distances, record scans of q‖qα‖‖qβ‖, the uniform
  quantity of section 8(a), and the p-adic quantity of section 8(c), in plain Python with no
  dependencies. The tests check Legendre's bounds and Borel's theorem on convergents, Lemma 3 on a
  dependent pair, the value at q = 41 for (√2, √3), and the cube root of 2 against OEIS A002945.

Nothing here claims progress on the conjecture itself.

## 12. Using the companion code

    python littlewood.py cf cbrt2 41                 # certified partial quotients of 2^(1/3)
    python littlewood.py records sqrt2 sqrt3 100000  # new minima of q ||q sqrt2|| ||q sqrt3||
    python littlewood.py value sqrt2 sqrt3 41        # certified interval for one value
    python littlewood.py uniform sqrt2 sqrt3 100000  # T * min_{q<=T} ||q x|| ||q y||, T = 1, 2, 4, ...
    python littlewood.py mixed sqrt2 2 100000        # new minima of q |q|_2 ||q sqrt2||
    python -m unittest test_littlewood               # or: python test_littlewood.py

Numbers: `sqrt2`, `sqrt(7)`, `cbrt2`, `2^1/3`, `2^2/3`, `golden`, and `surd:a,b,d,c` for
(a + b√d)/c. Python 3.6+, no dependencies.

## Sources

- J. W. S. Cassels, H. P. F. Swinnerton-Dyer, On the product of three homogeneous linear forms and indefinite ternary quadratic forms, Phil. Trans. Roy. Soc. A 248 (1955) 73–96. doi:10.1098/rsta.1955.0010
- L. G. Peck, Simultaneous rational approximations to algebraic numbers, Bull. AMS 67 (1961) 197–201.
- P. X. Gallagher, Metric simultaneous diophantine approximation, J. London Math. Soc. 37 (1962) 387–390. doi:10.1112/jlms/s1-37.1.387
- A. D. Pollington, S. L. Velani, On a problem in simultaneous diophantine approximation: Littlewood's conjecture, Acta Math. 185 (2000) 287–306. doi:10.1007/BF02392812
- G. A. Margulis, Problems and conjectures in rigidity theory, in Mathematics: Frontiers and Perspectives, AMS, 2000, 161–174.
- B. de Mathan, Conjecture de Littlewood et récurrences linéaires, J. Théor. Nombres Bordeaux 15 (2003) 249–266. doi:10.5802/jtnb.401
- B. de Mathan, O. Teulié, Problèmes diophantiens simultanés, Monatsh. Math. 143 (2004) 229–245. doi:10.1007/s00605-003-0199-y
- B. Adamczewski, Y. Bugeaud, On the Littlewood conjecture in simultaneous Diophantine approximation, J. London Math. Soc. 73 (2006) 355–366. doi:10.1112/S0024610706022617
- M. Einsiedler, A. Katok, E. Lindenstrauss, Invariant measures and the set of exceptions to Littlewood's conjecture, Ann. of Math. 164 (2006) 513–560. doi:10.4007/annals.2006.164.513 (arXiv:math/0612721)
- M. Einsiedler, D. Kleinbock, Measure rigidity and p-adic Littlewood-type problems, Compositio Math. 143 (2007) 689–702. doi:10.1112/S0010437X07002801
- A. Venkatesh, The work of Einsiedler, Katok and Lindenstrauss on the Littlewood conjecture, Bull. AMS 45 (2008) 117–134. doi:10.1090/S0273-0979-07-01194-9
- T. Gowers, Problems related to Littlewood's conjecture, blog post, 2009-11-17: https://gowers.wordpress.com/2009/11/17/problems-related-to-littlewoods-conjecture-2/
- F. Maucourant, A non-homogeneous orbit closure of a diagonal subgroup, Ann. of Math. 171 (2010) 557–570. doi:10.4007/annals.2010.171.557
- U. Shapira, A solution to a problem of Cassels and Diophantine properties of cubic numbers, Ann. of Math. 173 (2011) 543–557. doi:10.4007/annals.2011.173.1.11
- Y. Bugeaud, N. Moshchevitin, Badly approximable numbers and Littlewood-type problems, Math. Proc. Cambridge Philos. Soc. 150 (2011) 215–226. doi:10.1017/S0305004110000605
- D. Badziahin, On multiplicatively badly approximable numbers, Mathematika 59 (2013) 31–55. doi:10.1112/S0025579312000095
- Y. Bugeaud, Around the Littlewood conjecture in Diophantine approximation, Publ. Math. Besançon (2014) 5–18. doi:10.5802/pmb.1
- D. Badziahin, Computation of the infimum in the Littlewood conjecture, Exp. Math. 25 (2016) 100–105. doi:10.1080/10586458.2015.1031356
- F. Adiceam, E. Nesharim, F. Lunnon, On the t-adic Littlewood conjecture, Duke Math. J. 170 (2021) 2371–2419. doi:10.1215/00127094-2020-0077 (arXiv:1806.04478)
- S. Chow, A. Zafeiropoulos, Fully inhomogeneous multiplicative Diophantine approximation of badly approximable numbers, Mathematika 67 (2021) 639–646. doi:10.1112/mtk.12095
- S. Chow, L. Yang, An effective Ratner equidistribution theorem for multiplicative Diophantine approximation on planar lines, Invent. Math. 235 (2024) 973–1007. doi:10.1007/s00222-023-01233-1
- S. Chow, N. Technau, Dispersion and Littlewood's conjecture, Adv. Math. 447 (2024) 109697. doi:10.1016/j.aim.2024.109697 (arXiv:2307.14871)
- C. Stefanescu, Adv. Math. 461 (2025) 110062. doi:10.1016/j.aim.2024.110062
- D. Badziahin, Estimating lower limit in the p-adic Littlewood conjecture, arXiv:2506.16860 (2025).
- Vitorino, Vucusic, Some bounds related to the 2-adic Littlewood conjecture, Mathematika (2026). doi:10.1112/mtk.70073 (arXiv:2506.04110)
- P. Bandi, R. Fregoli, D. Kleinbock, Submanifold-genericity of ℝ^d-actions and uniform multiplicative Diophantine approximation, arXiv:2504.02258 (2025).
- R. Fregoli, D. Kleinbock, J. Number Theory 278 (2026) 570–621. doi:10.1016/j.jnt.2025.05.001
- S. Garrett, S. Robertson, Counterexamples to the p(t)-adic Littlewood conjecture over small finite fields, Math. Comp. 95 (2026) 1961–1986. doi:10.1090/mcom/4104 (arXiv:2405.14454)
- S. Robertson, Mathematika 72 (2026) e70064. doi:10.1112/mtk.70064
- T. Matala-Aho, T. Törmä, M. Vapa, Numerical upper bounds in the Littlewood conjecture, Exp. Math. (online 2026-04-02). doi:10.1080/10586458.2026.2615945
- J. Schleischitz, Disproof of the uniform Littlewood conjecture, arXiv:2603.12611 (2026).
- N. Moshchevitin, A note on uniform version of Littlewood inequality and Fibonacci numbers, arXiv:2605.26188 (2026).
- L. Lai, J. Sprang, On the P(t)-adic Littlewood conjecture in odd characteristics, arXiv:2606.00633 (2026).
- D. Badziahin, V. Pavlenkov, E. Zorin, Positive logarithmic Hausdorff measures of exceptional sets for the p-adic and t-adic Littlewood conjectures, arXiv:2608.22078 (2026).
- N. Shulga, The uniform Littlewood conjecture fails on a set of positive Hausdorff dimension, arXiv:2608.25059 (2026).
- V. Neckrasov, C. Wu, B. Yang, Winning property of counterexamples to Uniform Littlewood's Conjecture, arXiv:2608.24401 (2026).
- K. Dhanda, J. Flynn, A. Haynes, A geometric proof of Peck's theorem, arXiv:2609.29469 (2026).
- Y. Bugeaud, B. de Mathan, Refinements of Peck's theorem on simultaneous approximation to algebraic numbers, arXiv:2609.29360 (2026).
- OEIS A002945, Continued fraction for cube root of 2: https://oeis.org/A002945
- Wikipedia, Littlewood conjecture (edited 2026-10-05): https://en.wikipedia.org/wiki/Littlewood_conjecture
- Justin Sun Prize, MathDB pilot register: https://github.com/TheJustinSunPrize/awards/blob/main/problems/mathdb-pilot.md

## Licence

Text: CC BY 4.0. Code: Apache-2.0 (see LICENSE). Cite as: MV2 (Munim, Inc.), "The Littlewood conjecture in
October 2026", 2026.
