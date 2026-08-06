# Smooth-iteration computation in Mandelbrot-Explorer

This document describes how the engine computes the **smooth (continuous) iteration
count** — the scalar value that every non-trap coloring mode is built on — and every
optimization layered on top of it: **perturbation**, **series approximation (SA)**,
**Zhuoran rebasing**, **BLA**, the **floatexp** deep-zoom kernel, and the derived
quantities (EDE / normal map / DE overlay / stripe average).

All formulas below are the ones actually implemented in
`src/engine/mandel_perturbation.cpp`; line references point at the source.

---

## 0. Notation and the base iteration

The Mandelbrot iteration for a point $c\in\mathbb{C}$ is

$$
z_0 = 0,\qquad z_{n+1} = z_n^2 + c .
$$

A point *escapes* at the first $n$ with $|z_n| > R$, the **escape radius**
$R = 10^8$ (`_ESCAPE_RADIUS`, `mandel_perturbation.h:176`). Points that never
escape within `mxit` iterations are classified *interior* (sentinel `-2`).

Symbols used throughout:

| symbol | meaning |
|---|---|
| $C$ | the **reference** point (a full-precision GMP orbit) |
| $Z_n$ | reference orbit, $Z_{n+1}=Z_n^2+C$ |
| $c=C+\delta c$ | a pixel, offset $\delta c$ from the reference |
| $z_n=Z_n+\delta z_n$ | the pixel orbit, tracked via its delta $\delta z_n$ |
| $J_n=\mathrm{d}z_n/\mathrm{d}c$ | total derivative (for EDE / normal / DE) |
| `dx`,`dy` | per-pixel step in the complex plane |

---

## 1. The smooth iteration count

### 1.1 Formula

On escape at iteration $n$ with squared modulus $r = |z_n|^2$ (`zrad`), the engine
returns (e.g. `mandel_perturbation.cpp:1523`, `:2052`, `:360`):

```c
mu = n + 1 - log( log(zrad)/2 / log(2) ) / log(2)
```

which is

$$
\boxed{\;\mu = n + 1 - \log_2\!\big(\log_2 |z_n|\big)\;}
$$

### 1.2 Derivation

Near escape the quadratic term dominates, $z_{n+1}\approx z_n^2$, so
$\log|z_{n+1}| \approx 2\log|z_n|$, i.e. the quantity

$$
s_n \;=\; \log_2\log|z_n|
$$

increases by **exactly 1** per iteration once $|z|$ is large. Therefore $n - s_n$
is (asymptotically) constant along the orbit, and the *fractional overshoot*
past the bailout is $s_n \bmod 1$. Defining the continuous count as
$\mu = n + 1 - s_n$ makes $\mu$ increase smoothly and continuously as the pixel
crosses a bailout contour, removing the integer "banding" of raw escape time.
Because $R=10^8$ is large, the additive constant that formally depends on $R$ is
absorbed into a palette phase offset and is dropped.

Rewriting the code: with $r=|z|^2$,

$$
\tfrac{1}{2}\ln r = \ln|z|,\qquad
\frac{\tfrac12\ln r}{\ln 2} = \log_2|z|,\qquad
\frac{\ln(\log_2|z|)}{\ln 2} = \log_2\log_2|z|,
$$

which is exactly the C expression. (Using $\log_2|z|$ rather than $\log|z|$ only
shifts $\mu$ by the constant $\log_2\ln 2$.)

Everything below exists to compute the escaping $n$ and $r$ (and, for the derived
modes, $J_n$) **fast and accurately at extreme zoom**.

### 1.3 Custom quadratic deep-zoom compatibility

Custom frames stay on the direct expression AVX2 path through zoom $10^{12}$.
Above that boundary they use this quadratic perturbation/SA/BLA engine only when
the compiled source and specialized runtime are the exact recognized `z*z+c`
recurrence in the c-plane, fixed `z0` is positive zero in both components, and the
Custom bailout is finite, at least 1, and has a finite squared value. A scoped
per-frame override makes the reference orbit, SA/BLA, perturbation loops, and
direct fallbacks use that bailout without changing ordinary Mandelbrot frames.
The whole frame must also have $z_1=c$ inside the bailout disk because the
optimized engine starts testing at $z_2$.

The two renderers historically expose slightly different numeric conventions.
The deep Custom writer therefore adds the constant
$-\log(\log 2)/\log 2$ for Smooth and halves Distance (the production engine uses
$\log |z|^2$ where Custom uses $\log |z|$). Feather and OrbitTrap are not routed
in stage 1 because Custom includes $z_1=c$ in their accumulators while the
production quadratic accumulators begin at $z_2$; Raw also remains direct. Those
modes and every incompatible formula retain the $10^{12}$ cap.

---

## 2. Perturbation

### 2.1 Why

At zoom $10^{40}$ the pixel spacing is $\sim 10^{-40}$; `double` (≈16 digits) cannot
even represent two adjacent pixels distinctly. Computing every pixel in GMP is
correct but ~100× too slow. **Perturbation** computes one expensive reference orbit
$Z_n$ in high precision, then each pixel as a small `double` **delta** relative to it.

### 2.2 The delta recurrence

Write $c = C + \delta c$ and $z_n = Z_n + \delta z_n$. Substituting into
$z_{n+1}=z_n^2+c$:

$$
Z_{n+1}+\delta z_{n+1} = (Z_n+\delta z_n)^2 + (C+\delta c)
= \underbrace{Z_n^2+C}_{Z_{n+1}} + 2Z_n\,\delta z_n + \delta z_n^2 + \delta c .
$$

Cancelling $Z_{n+1}$ gives the fundamental **perturbation recurrence**

$$
\boxed{\;\delta z_{n+1} = 2\,Z_n\,\delta z_n + \delta z_n^{\,2} + \delta c\;}
$$

implemented at `mandel_perturbation.cpp:1485–1495`:

```c
// dz = dz^2 + 2*dz*z + dc
tmp   = (dzr*_zfr[rkm1] - dzi*_zfi[rkm1])*2 + dzr2 - dzi2 + dcr;
dzi   = (dzr*_zfi[rkm1] + dzi*_zfr[rkm1])*2 + dzr*dzi*2 + dci;
dzr   = tmp;
```

The reconstructed orbit and its modulus are

$$
z_n = Z_n + \delta z_n,\qquad r = |z_n|^2
$$

(`:1497–1499`). The key point: $\delta z_n$ and $\delta c$ are $O(\text{pixel size})$
in *magnitude* but each individually well-scaled, so `double` arithmetic on the
delta is accurate — even though $z_n$ itself is far below `double`’s spacing.

The reference $Z_n$ is stored as a low-precision shadow `_zfr[k]/_zfi[k]`
(`long double` cast, `:2388–2390`); only the *pixel* uses `double`.

### 2.3 The reference orbit (GMP), cheaply

`calCoefficient` (`:2368`) advances the reference in GMP with two tricks:

* **Rotating two-buffer orbit.** `_z_re[i&1]` / `_z_re[(i-1)&1]` alternate parity so
  current and previous are always distinct without copying (`:2376`).
* **Karatsuba complex square.** $z^2=(a+bi)^2$ needs only **two** big multiplies
  instead of the naive three:

  $$
  \Re = (a+b)(a-b),\qquad \Im = 2ab .
  $$

  (`:2379–2385`.) The reference build is the only inherently serial per-frame cost
  ($\sim 0.4$ s at $10^{432}$), so its multiply count matters.

Precision is set to $|\text{exp}_2(\text{scale})| + 30$ bits, just enough to resolve
the view.

---

## 3. Series Approximation (SA)

### 3.1 Idea

For the *first* many iterations, $\delta z_n$ is an analytic function of $\delta c$
and can be written as a power series. Evaluating that series lets us **skip** those
iterations entirely (per pixel) and jump straight to iteration `SA_it`.

Let $\delta = $ `_SA_delta` (the corner delta, a fixed $O(\text{half-view})$ scale,
`:2326–2331`) and $u = \delta c/\delta$. The engine represents

$$
\delta z_n \;=\; \sum_{k=0}^{N-1} A_k^{(n)}\; u^{\,k+1}.
$$

### 3.2 Coefficient recurrence

Substitute the series into $\delta z' = 2Z\,\delta z + \delta z^2 + \delta c$ and
equate powers of $u$. The $\delta z^2$ term is a Cauchy (convolution) product, so

$$
A_j' \;=\; 2 Z\,A_j \;+\!\!\sum_{k+l=j-1}\!\! A_k A_l \;+\; [\,j=0\,]\,\delta .
$$

The code halves the convolution using its symmetry
($\sum_{k<j/2}2A_kA_{j-1-k}$ plus a middle square when $j$ is odd),
`:2418–2435`:

```c
_Adf_new[j] = 2 * z[i-1] * _Adf_old[j];
for (k = 0; k < j/2; ++k) _Adf_new[j] += 2 * A_old[k] * A_old[j-k-1];
if (j % 2)               _Adf_new[j] +=     A_old[j/2] * A_old[j/2];
...
_Adf_new[0] += _SA_delta;                 // the +delta seed of the linear term
```

### 3.3 Validity test

The truncated series is trustworthy only while the dropped tail is negligible.
`SACheckMagnitude` (`:2455`) computes, over the $N$ coefficient magnitudes
(as base-2 exponents), the running **prefix-min** and **suffix-max**, and accepts
truncation order $i$ at the first place where

$$
\min_{k\le i}\exp(A_k)\;-\;\max_{k> i}\exp(A_k)\;\ge\;80\ \text{bits},
$$

i.e. the high-order terms are $\ge 2^{80}$ smaller than the retained ones
(`:2476–2478`). If no such gap exists the series has stopped converging and SA is
disabled (`_SA_flag=false`). The last valid iteration is `_SA_it`, the truncation
order `_SA_order`.

### 3.4 Seeding the delta loop

Each pixel evaluates the series by Horner’s scheme (`:1355–1359`, `:1400–1404`):

```c
dz = 0;
for (x = _SA_order; x >= 0; --x) { dz += _Adf_old[x]; dz *= dc; dz /= _SA_delta; }
```

giving $\delta z_{\text{SA\_it}} = \sum_k A_k (\delta c/\delta)^{k+1}$, and the
iteration then resumes at $j=\texttt{SA\_it}+1$. Under EDE a **second** series
$B_j$ (the derivative $\mathrm{d}\,\delta z/\mathrm{d}c$) is carried the same way
(`:2426–2434`).

---

## 4. Zhuoran rebasing (glitch elimination)

### 4.1 The glitch

Perturbation reconstructs $z_n = Z_n + \delta z_n$. If the reference orbit passes
**close to zero** ($|Z_n|$ tiny) while the pixel’s $|\delta z_n|$ is not, then
$Z_n \approx -\delta z_n$ and the reconstruction suffers catastrophic cancellation —
the linear driver $2Z_n\delta z_n$ loses its significant digits and the pixel
“glitches’’ (wrong color blobs).

### 4.2 The fix

**Zhuoran’s rebasing**: whenever the true orbit value is *smaller* than the delta,
the delta is no longer the small, well-conditioned quantity — $z_n$ itself is. So
restart the perturbation from the reference’s **start** (index $0$, the critical
point) using the reconstructed $z_n$ as the new delta. Condition and action
(`:1561–1574`):

$$
|z_n|^2 < |\delta z_n|^2 \;\;\Longrightarrow\;\;
\delta z_n \leftarrow z_n,\quad k \leftarrow 0 .
$$

```c
if (zrad < dzr*dzr + dzi*dzi || (!per && k == mx_ref_it)) {
    ...
    dzr = zr; dzi = zi;          // dz <- reconstructed z
    k = 0;                       // reference index back to the start
}
```

Because the *new* $|\delta z| = |z_n|$ is by construction the smaller of the two,
the delta stays the minimal representation at every step, and the cancellation
never happens. This single rule makes the double path essentially glitch-free — on
the floatexp path it is enough on its own (`:1320–1324`, every pixel resolves in one
pass).

### 4.3 Two more rebase triggers

* **Reference exhausted.** A non-periodic reference of length `mx_ref_it` also
  rebases at $k=\texttt{mx\_ref\_it}$ (wrap for periodic references), `:1561`.
* **Hard glitch guard (Pauldelbrot).** If a rebase would incur a large magnitude
  jump, $|\delta z|^2/|z|^2 > 10^7$, precision is already lost; the pixel is flagged
  and a **new reference** is requested (`:1561–1566`). The next reference is chosen
  glitch-guided as the flagged pixel with the smallest orbit magnitude
  $|z|^2/|Z_{\text{ref}}|^2$ (`:1588–1595`, `createRef` `:2280`).

---

## 5. BLA — Bilinear (linear) Approximation

### 5.1 Idea

When $|\delta z|$ is small, the quadratic term $\delta z^2$ in
$\delta z' = 2Z\delta z + \delta z^2 + \delta c$ is negligible and one step is
**affine** in $(\delta z,\delta c)$:

$$
\delta z_{n+1} \approx 2 Z_n\,\delta z_n + \delta c .
$$

A run of $L$ such steps composes into a single affine map

$$
\boxed{\;\delta z_{n+L} = A\,\delta z_n + B\,\delta c\;}
$$

so we can **skip $L$ reference iterations at once**.

### 5.2 Coefficients and merging

Single step (level 0) at reference index $s$: $A = 2Z_s$, $B = 1$ (`:1640–1647`).
Composing map $x$ (skip $l_x$) followed by map $y$:

$$
A_z = A_y A_x,\qquad B_z = A_y B_x + B_y,\qquad l_z = l_x + l_y
$$

(`:1665–1668`). The engine builds a **binary tree**: level $p$ merges neighbour
pairs, so its entries skip $2^p$ steps (`:1654–1691`).

### 5.3 Validity radius

The skip is valid only where the dropped $\delta z^2$ is below rounding. Using
$\varepsilon = 2^{-53}$ and the worst-case pixel offset `dcmax` (farthest image
corner from the reference, `:1628–1631`), the single-step radius is the
mathr/Zhuoran bound (`:1642–1645`):

$$
R = \max\!\Big(0,\ \varepsilon\,\frac{|Z_s| - \texttt{dcmax}}{|2Z_s| + 1}\Big),
\qquad r^2 = R^2 .
$$

Merged levels take $R = \min\!\big(R_x,\ (R_y - |B_x|\,\texttt{dcmax})/|A_x|\big)$
(`:1680–1685`) so the whole run stays valid.

> **`dcmax` must use the actual reference position.** After glitch
> re-referencing the reference is off-centre, so a centre-based bound underestimates
> edge-pixel $|\delta c|$ and a skip can jump over an escape (this was the
> `flake@1e157` corner misclassification, `:1622–1627`).

### 5.4 Applying a skip

`tryBLA` (`:1700`) picks, for the current reference index $s$, the **largest aligned**
level whose radius contains the current $|\delta z|^2$. Levels start at index
$1 + i\,2^p$, so $s$ is aligned at levels $0..\text{tz}$ where $\text{tz}$ =
trailing zeros of $s-1$ — the search starts there, no per-level alignment test
(`:1705–1713`). It also **never skips past an escape**: it reconstructs
$z_{\text{land}} = Z_{\text{land}} + (A\delta z + B\delta c)$ and rejects the skip if
that already exceeds the bailout (`:1720–1726`).

BLA is disabled for orbit traps (which need every orbit point) and combined with
the derivative map $C,E$ under EDE (`:1728–1734`). At $10^{432}$ BLA skips ~99.5 %
of all reference iterations (avg skip ≈ 1220), turning the deep delta loop from
$\sim$20 s into $\sim$2 s per frame.

---

## 6. The floatexp deep-zoom kernel

### 6.1 Why

Past scale $\sim 10^{300}$ the pixel delta $\delta z$ (and intermediates like
$2Z\delta z$) **underflow `double`** (min normal $\approx 10^{-308}$). The engine
escalates to a **rescaled floatexp** path — a hard gate at `scale > 1e280`, plus a
*content-aware* escalation whenever the just-built reference passes close enough to
zero that $2Z\,\delta z$ would denormalize (`:878–891`).

`FloatExp` (`src/engine/floatexp.h:17`) is a `double` mantissa $m$
($\tfrac12\le|m|<1$) with a 64-bit base-2 exponent $e$ — value $m\cdot 2^{e}$, so its
exponent range is effectively unbounded while arithmetic stays as cheap as a couple
of `double` ops.

### 6.2 Rescaled delta

Carry the delta as $\delta z = S\cdot w$, where $S$ is a **FloatExp scale** and
$w=(w_r,w_i)$ is a **`double`** unit-ish vector; likewise $d = \delta c/S$ and the
scalar $s = |S|$ in `double`. Substituting $\delta z = Sw$ into the perturbation
recurrence and dividing by $S$ gives the rescaled step (`pixelRescaled` `:2001–2007`):

$$
w' = 2 X_m\,w + s\,w^2 + d,\qquad X_m = Z_m\ (\text{double shadow}),
$$

with reconstruction $z = X_m + S\,w$ (`:2025`). To keep the mantissa healthy, $w$ is
**periodically renormalized**: when $|w|^2$ leaves $[10^{-16},10^{16}]$, fold its
magnitude into $S$ ($S\leftarrow S|w|$, $w\leftarrow w/|w|$), `:2076–2084`. All the
*decisions* that need true magnitude — the BLA radius test (`tryBLAfe` `:1746`), the
Zhuoran rebase ($|z|^2<|\delta z|^2$), and the reference-end wrap — are done in
FloatExp (`:2086–2111`), while the hot quadratic step stays in `double`. This is
byte-identical whether run scalar or 4-wide SIMD; the deep path is
memory/latency-bound on the reference gather, so SIMD only ≈ matches scalar.

---

## 7. Derived quantities (all built on the same orbit)

Every mode below reuses the smooth engine and adds one accumulator.

### 7.1 Total derivative $J = \mathrm{d}z/\mathrm{d}c$

Differentiating $z_{n+1}=z_n^2+c$ w.r.t. $c$:

$$
\boxed{\;J_{n+1} = 2 z_n J_n + 1,\qquad J_0 = 0\;}
$$

(`pixelRescaled` `:2010–2013`, in its own rescaled floatexp $J=S_J\,(j_r,j_i)$ so it
survives $|J|\sim 1/\text{dx}$ overflow). On the double path the derivative is tracked
as a **perturbation** $dd = J - D$ of the reference derivative $D$
($D_{n+1}=2D_nZ_n+1$, `:2396–2412`):

$$
dd' = 2\big(dd\cdot z + \delta z\cdot(D+dd)\big)
$$

(`:1472–1483`), and it is carried through a BLA skip by the same linear map,
$J \to A J + B$ (`:1728–1734`, `:1983–1991`). On rebase, $dd$ resets consistently
with $\delta z$ (`:1569–1572`).

### 7.2 Exterior Distance Estimate (EDE)

The Milnor distance from $c$ to the set boundary is
$\mathrm{dist} = |z_n|\ln|z_n| / |z_n'|$ with $z_n' = J_n$. The engine returns it in
**pixel units** (÷`dx`), `:1509`, `:2030–2034`:

$$
\mathrm{dist}_{\text{px}} = \frac{|z_n|\,\ln|z_n|^2}{\text{dx}\,\cdot|J_n|}
= \frac{\sqrt{r}\,\ln r}{\text{dx}\,\sqrt{J_r^2+J_i^2}} .
$$

### 7.3 Normal map / relief

The surface normal angle for Lambert (slope) lighting is

$$
\theta = \arg(z_n) - \arg(J_n) = \arg\!\big(z_n/(\mathrm{d}z_n/\mathrm{d}c)\big),
$$

`:1513`, `:2038`. The **DE overlay** mode uses the same $J$ to draw the pixel-
normalized distance as a black-and-white filament layer over the smooth base
(`:1517`, `:2045`).

### 7.4 Stripe Average Coloring (SAC / "Feather")

Average the stripe function $t_n = \tfrac12 + \tfrac12\sin(7\,\arg z_n)$ along the
orbit, blended by the smooth fraction for continuity (`SacAccum` `:48–78`). Two
subtleties make it deep-zoom-clean:

* **Windowed tail (Bartlett).** Only the last $W=256$ iterations are averaged, with
  **triangular** weights (newest weight $W$, fading linearly to 0 at age $W$). The
  finite support keeps it reference-clean; the linear fade avoids a hard far-edge
  cliff (`:40–41`, `:58–76`). The full-orbit average goes flat past $\sim10^{250}$.
* **Skip handling.** A BLA skip breaks the tail window, so it is **reset**
  (`:1459`, `:1978–1979`); for the classic full-average variant the omitted stripe
  contributions are restored from a reference-orbit **prefix sum** `_sacRefPre`
  ($z\approx X$ during a valid skip), avoiding a pan-dependent radial halo
  (`:905–916`, `:1980–1981`).

The value blends the windowed average $a_1$ and the average-without-newest $a_2$ by
the smooth fraction (`:67–76`):
$\text{value} = a_2 + (a_1-a_2)\,\mathrm{frac}$.

### 7.5 Orbit trap

Track the orbit’s closest approach to a composite trap — a point at $0$, the
Pickover cross ($x{=}0$ or $y{=}0$), and the circle $|z|{=}0.5$ — and map
$-\log_{10} d$ (clamped $\ge 0$) plus a small smooth-count term to the palette
(`TrapAccum` `:83–108`). No angle term is added: a wrapping $\arg$ on a linear
palette coordinate seams at the $\pm\pi$ branch cut (`:102–105`). BLA is disabled
here (every orbit point is needed).

---

## 8. Other optimizations (orchestration)

`Mandel::Compute` (`:663`) wires the above together:

* **Shallow vs deep.** `scale > 1e6` switches from the direct `double` iteration
  (`floatPointCompute` `:307`, also the ground-truth/oracle path) to
  perturbation+SA+rebase (`:708–711`).
* **Periodic reference.** If a minibrot **nucleus** of period $p$ sits in view
  (`findNucleus`), build **one period** in GMP and index it modulo $p$ in the delta
  loop, instead of a full `mxit`-long orbit — a big serial-cost win, auto-enabled on
  the floatexp deep path (`:820–862`). Requires BLA on, non-EDE, non-SAC.
* **Interior detection.** Brent-style periodicity (a "tortoise" saved at geometric
  intervals) proposes a candidate cycle; it is **confirmed** over several periods by
  (a) close returns and (b) multiplier $\prod 4|z|^2 < 1$ (attracting) → interior,
  short-circuiting the walk to `mxit` (`:1528–1558`, `:2056–2074`). Auto-gated off if
  a coarse probe finds no interior pixel (`:936–942`). A BLA skip invalidates the
  detector state, so it restarts after a skip (`:1460–1467`).
* **Coarse-to-fine preview.** A strided ($1/16$-work) pass paints a blocky preview
  picked up immediately by the async display, then the full pass sharpens in place
  (`:728–747`, `:920–962`).
* **SIMD (AVX2, 4-wide).** `solveShallowSimdList`, `solveSimd4` (`:1790`) and
  `solveRescaledSimd4` (`:2121`) process 4 pixels/lanes with the reference read via
  per-lane gather; escape/rebase/glitch are a short scalar pass over the 4 lanes.
  Written in the same op order (no FMA contraction) so output is **bit-identical** to
  scalar. SIMD is used only where the quadratic step is the whole cost (BLA/EDE/SAC/
  trap/interior paths run scalar, `:1328`, `:1263`).
* **Multi-pass re-referencing.** Unresolved (glitched) pixels remain in the pending
  set `s`; each pass picks a new glitch-guided reference, rebuilds BLA, and re-runs
  `stepParallel` until `s` is empty (`:975–999`).

---

## 9. Correctness harness

Any change to the engine is validated against a GMP oracle:

```
build\verify.exe flake 128 96      ->  checksum(pert)=0x57b9524a, class mismatch: 0, PASS
```

`accuratePointCompute` (`:457`) is the per-pixel GMP ground truth;
`render.exe MANDEL_ORACLE` renders whole frames in GMP for image diffs. The
perturbation/SA/rebase/BLA/floatexp stack is required to match the oracle to within
smooth-count rounding, with **zero** interior/exterior class mismatches.

---

### Reference-line index (quick map)

| topic | function / lines |
|---|---|
| smooth formula | `:1523`, `:2052`, `:360` |
| perturbation step | `stepParallel` `:1485–1499` |
| GMP reference (Karatsuba) | `calCoefficient` `:2379–2385` |
| SA coefficients | `calCoefficient` `:2418–2435`; check `SACheckMagnitude` `:2455` |
| SA seeding | `:1355–1359`, `:1400–1404` |
| Zhuoran rebase | `:1561–1574`; floatexp `:2086–2111` |
| BLA build | `buildBLA` `:1616–1695` |
| BLA apply | `tryBLA` `:1700`; floatexp `tryBLAfe` `:1746` |
| floatexp kernel | `pixelRescaled` `:1909` |
| total derivative $J$ | `:2010–2013`; delta form `:1472–1483` |
| EDE / normal / DE | `:1508–1519`, `:2030–2047` |
| SAC / Feather | `SacAccum` `:48–78` |
| orbit trap | `TrapAccum` `:83–108` |
| orchestration | `Compute` `:663` |
