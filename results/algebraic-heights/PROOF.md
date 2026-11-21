# Proof

## Prime logarithms at algebraic heights

Write $`K=\overline{\mathbf Q}\subset\mathbf C`$, $`a=\log 2`$ and $`K_0=\mathbf Q(\zeta_{16})^+`$. For $`i=(p,q)`$ equal to $`A=(547,727)`$ or $`B=(571,599)`$, the base is $`F_i=K_0(\sqrt{pq})`$ and the quadratic extension is $`E_i=K_0(\sqrt{-p},\sqrt{-q})`$. Put $`\epsilon_v=+1`$ for a split prime ideal of $`F_i`$ and $`-1`$ for an inert one. Omit primes above 2 and the host’s two conductors. The published host conventions also omit 343057911569 and 343070636977 on A, and 413718278401 and 603339156001 on B. Enlarging these omissions by a fixed finite set preserves the eventual result. [Trace on A](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [Trace on B](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean) · [Prime-ideal decomposition](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/PROOF.md#nine-frequencies-on-one-height-axis)

For $`b\gt1`$ put $`B_b=2\sqrt b(\sqrt b-1)^2`$ and $`W_b(t)=\max(0,\min(b(1-1/t),b^2/t-1))`$. The score is

```math
S_{i,b}(s)=-\frac1{B_b\sqrt s}\sum_v\epsilon_v\log(Nv)\,W_b(Nv/s).
```

Only norms between $`s`$ and $`b^2s`$ contribute. At positive algebraic $`s,b`$ with $`b\gt1`$, the coefficient of each rational prime logarithm is algebraic. Indeed $`Nv=\ell^f`$ and $`\log(Nv)=f\log\ell`$. Every such $`\ell`$ is odd. The real logarithms of any finite set of distinct rational primes, including 2, are linearly independent over $`\mathbf Q`$ by unique factorization after clearing denominators and exponentiating. Baker’s theorem promotes this to linear independence over $`K`$. Consequently the $`K`$-span of the normalized readings $`S_{i,b}(s)/a`$ at all algebraic inputs meets $`K`$ only at zero. Each reading is either zero or transcendental. [Baker’s theorem](https://doi.org/10.1112/S0025579300003971) · [Score conventions](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/README.md)

This also gives an exact test at a small algebraic start. Make a column for each score, with a row for each contributing odd rational prime. Enter the coefficient of $`\log\ell`$ in the unnormalized score, adding all ideals above $`\ell`$. For a relation among 1 and the normalized readings, multiplying by $`a`$ forces the constant coefficient to vanish separately. The other coefficients of the relation belong to the nullspace of this finite matrix, by Baker’s theorem. Conversely any vector in that nullspace gives an identity.

## The rank at algebraic starts

Fix a finite list of observations $`(i,b,u)`$, with real algebraic $`b\gt1`$ and $`u\gt0`$. Widths, scales and observations may repeat. Take $`S_{i,b}(uy)/a`$ on the list at a common positive real start $`y`$. In each host $`i`$ form a matrix $`J_i`$ with one row per distinct knot. Put the following entries in each column and fill the other rows with zero.

| Knot | Entry multiplied by $`B_b\sqrt u`$ |
| --- | --- |
| $`u`$ | $`b`$ |
| $`bu`$ | $`-(b+1)`$ |
| $`b^2u`$ | 1 |

**Theorem.** For every sufficiently large positive algebraic $`y`$ the linear relations over $`K`$ among 1 and the readings have zero constant coefficient and their remaining coordinates form exactly $`\ker J_A\oplus\ker J_B`$. The dimension of the span, including 1, is $`1+\mathrm{rank}(J_A)+\mathrm{rank}(J_B)`$. Every relation from the two kernels holds at *every* real $`y\gt0`$. The onset depends on the list. It works for every choice of coefficients for that list.

Let $`z`$ denote norm divided by $`y`$. After removing $`1/\sqrt y`$, the score kernel is $`W_b(z/u)/(B_b\sqrt u)`$. Multiplication by the common positive variable $`z`$ gives the piecewise-linear numerator

```math
T_{b,u}(z)=\max(0,\min(b(z-u),b^2u-z)).
```

It is zero before $`u`$ and after $`b^2u`$. Its slopes between the three knots are $`b`$ and $`-1`$. A combination of the scaled hats vanishes identically exactly when its slope jumps all vanish, since it starts at zero with slope zero. This proves that every vector in the kernel gives an identity at all real starts. It works at each actual norm, including prime powers.

For the reverse inclusion, select two reduced classes modulo $`M=16\cdot547\cdot727\cdot571\cdot599`$. For A prescribe 1 modulo 16, 547, 727 and 599, and a quadratic nonresidue modulo 571. For B prescribe 1 modulo 16, 727, 571 and 599, and a nonresidue modulo 547. The Chinese remainder theorem supplies the classes. A prime from the class for $`i`$ splits in $`K_0`$ and $`F_i`$, and all eight degree-one ideals have the same sign in $`E_i/F_i`$. In the other base field the residue degree is two. In fact the degree-one criterion is splitting in $`K_0`$ together with equal quadratic symbols at the two conductors. The signs follow by quadratic reciprocity. [Local splitting formula](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/PROOF.md#nine-frequencies-on-one-height-axis) · [Reduced-class checks](verify_rank.py)

Choose two distinct primes from the class for $`i`$ with $`\ell/y`$ in each open interval between successive knots of $`J_i`$. For every fixed $`0\lt c\lt d`$ the number in $`(cy,dy)`$ is asymptotic to $`(d-c)y/(\varphi(M)\log y)`$ by the prime number theorem in arithmetic progressions. There are finitely many intervals and classes. Thus this choice is possible for *every* large enough real $`y`$. All selected primes can avoid the omissions, with $`\ell^2`$ exceeding the largest support endpoint in the entire list. The latter follows because the smallest knot is a positive constant and the largest norm is a fixed multiple of $`y`$. An empty list for one host needs no primes there. [Prime number theorem in fixed progressions](https://kskedlaya.org/ant/chapter-5.html)

Now suppose a numerical relation holds at such an algebraic $`y`$. Each prime log has zero coefficient by Baker’s theorem. At a selected prime, only the degree-one row in the chosen host contributes. Its signed multiplicity and all common normalizations are nonzero. Dividing them out shows that the combined hat in this host vanishes at both chosen points of each interval. An affine function with two distinct zeros vanishes on the entire interval. Continuity supplies the knot values. Hence the relation lies in that host’s kernel. The two arguments together prove the theorem.

With a fixed width per host, distinct dilations have independent columns. In a putative dependence take the smallest left endpoint among the columns with nonzero coefficients. Nothing else can cancel its first slope jump. Thus $`m`$ distinct common dilations give rank $`2m+1`$ including 1. For $`u=1,\sqrt2,2`$ this also gives independence of 1, the two midpoint readings $`S_{i,2}(\sqrt2 y)/a`$ and both moves $`D_i(y)=(S_{i,2}(2y)-S_{i,2}(y))/a`$. The four observable rows are independent combinations of the six score rows.

## A finite range for three observations

For $`u=1,\sqrt2,2`$ and width two, use three prime bands $`(1.1y,1.4y)`$, $`(1.5y,1.9y)`$ and $`(2.1y,2.6y)`$. Choose one prime per band and field with $`\ell\equiv\pm1\bmod16`$, equal symbols in that field and unequal symbols in the other. Either common sign works. The three rows are lower triangular in the columns ordered by dilation. Each selected prime is in the positive support of its diagonal hat and below every later left endpoint. All squares exceed $`8y`$ on the checked range. The six rows therefore give two invertible blocks.

A prime $`p`$ in a band $`(Ly/10,Ry/10)`$ covers the open interval $`(10p/R,10p/L)`$ of starts. Six chains, with lengths 28, 29, 32, 29, 30 and 32, overlap strictly and cover both endpoints of $`[1500,1000000]`$. The checker first sieves and uses Euler’s criterion, then rechecks every selected prime by trial division and every sign from explicit square-residue tables. Integer comparisons prove the overlaps. Thus the full rank holds at every algebraic height in this closed range. [Finite arithmetic](verify_rank.py) · [Transcript](evidence/check.txt)

## An explicit infinite range and a prime bound

The six readings at scales $`1,\sqrt2,2`$ have full rank at **every algebraic $`y\ge e^{366}`$**. The proof gives a bound on ordinary primes as well. For a rational prime $`p\equiv\pm1\bmod16`$, let $`e_A(p)`$ and $`e_B(p)`$ be the two signed degree-one splitting coefficients divided by 8, with the split sign positive. A degree exceeding one gives zero. For a unit vector $`v\in\mathbf R^2`$ and $`L\ge e^{366}`$, let $`\mathcal P_\sigma(L)`$ contain the primes $`L\lt p\le4L/3`$ with $`p\equiv\pm1\bmod16`$ and $`\sigma v\cdot(e_A(p),e_B(p))\gt1/32`$. The prime bound is

```math
\sum_{p\in\mathcal P_\sigma(L)}\log p\gt L/864,\qquad \sigma=\pm1.
```

Each set therefore has more than $`L/(864\log(4L/3))`$ distinct primes. All the fixed conductor and host omissions lie below this range. Reversing both sign conventions leaves the two bounds intact.

### The real characters needed for an explicit onset

Write $`J`$ for the indicator of $`n\equiv\pm1\bmod16`$ and $`\chi_r(n)=(n/r)`$ for the quadratic symbol at a conductor $`r`$. At a good prime in these classes the two components are $`(\chi_{547}+\chi_{727})/2`$ and $`(\chi_{571}+\chi_{599})/2`$. Their squares and product use only pairs of the four characters. The required moduli are 16, $`16r`$ and $`16rs`$ for distinct conductor primes. The largest is 6,967,568.

Every primitive real character inducing a nonprincipal real character at one of these moduli has fundamental discriminant $`D=d\prod_{r\in I}(-r)`$ with $`d\in\{1,-4,8,-8\}`$ and $`I`$ a subset of at most two of the four primes. Excluding the principal character leaves 43. This follows from the Chinese remainder theorem, the unique quadratic character at each odd prime, quadratic reciprocity, and the fact that the real characters modulo 16 all factor through 8. An imprimitive Euler factor at positive real $`s`$ cannot vanish. The principal factor has a pole at 1 and has no real zero in $`(0,1)`$ by the alternating zeta series.

The certificate gives $`L(s,\chi_D)\gt0.363`$ for **all 43** on the entire real segment $`15/16\le s\le1`$. Write $`A(n)=\sum_{j\le n}\chi_D(j)`$ and bound $`|A(n)|`$ by the exact maximum over a full period. Periodicity and the vanishing period sum justify

```math
L(s,\chi_D)=s\int_1^\infty A(\lfloor t\rfloor)t^{-s-1}\,dt.
```

At $`N=65536`$ the tail is at least $`-\max|A|\,N^{-s}`$. Divide $`1-s`$ into the sixteen cells $`[k/256,(k+1)/256]`$. On each of 1,002 integer-ended integration blocks $`[l,h)`$ take the minimum partial sum $`m`$. If $`m\ge0`$, a lower bound for the block is $`m(255-k)l^{k/256}(1/l-1/h)/256`$. For $`m\lt0`$ use $`m(256-k)h^{(k+1)/256}(1/l-1/h)/256`$. These follow by writing $`t^{-s-1}=t^{-2}t^{1-s}`$. Eight directed integer square roots bound each 256th root. Integer floors in every block and the upward-rounded tail produce a lower bound for the full integral. [Real-axis certificate](real_characters.py)

Bordignon’s Theorems 1.1 and 1.2 isolate a possible real exceptional zero at $`\beta\ge1-1/(6.3970\log q)`$. This exceeds $`15/16`$ for $`q\ge16`$, so the exceptional term vanishes for every modulus used here. In Table 6 of the **published version**, the row $`(Y_0,\alpha_1,\alpha_2,C)=(5.9,3,1,7.58)`$ applies at $`X\ge e^{366}`$ since $`q\le366^3`$ and $`e^{5.9}\lt366`$. The row $`(5.8,2,2,0.27)`$ also applies to 16 and each $`16r`$, all below $`366^2`$. In a reduced class the error about $`X/\varphi(q)`$ for the von Mangoldt sum is bounded by $`CX/(\varphi(q)(\log X)^{\alpha_2})`$. [Bordignon, Theorems 1.1–1.2](https://nyjm.albany.edu/j/2021/27-54v.pdf#page=3) · [Published Table 6](https://nyjm.albany.edu/j/2021/27-54v.pdf#page=21)

### Two moments of the local traces

Exactly one quarter of the reduced classes have $`J=1`$. On these classes the average of any single $`\chi_r`$ or pair $`\chi_r\chi_s`$ used below vanishes by the Chinese remainder theorem. Summing the errors over the selected classes costs a quarter of the corresponding full bound. Removing proper prime powers costs at most $`\sqrt X(\log X)^2/\log2`$, by counting at most $`\sqrt X`$ possible bases for each exponent. This is below $`0.005 X/(\log X)^2`$ for $`X\ge e^{366}`$. Indeed $`t^4e^{-t/2}`$ decreases for $`t\ge366`$ and $`183^{10}/10!\gt400\cdot366^4`$, with $`\log2\gt1/2`$. Thus each pair moment has absolute endpoint error at most $`1.90X/\log X`$. The constant and single-character moments have error at most $`0.073X/(\log X)^2`$. The constant moment has main term $`X/4`$, the signed moments have main term zero. Each pair uses its own modulus, and any prime dividing the relevant odd modulus has zero character value.

Take $`L\ge e^{366}`$, $`U\gt L`$, $`\Delta=U-L`$ and $`U+L\le7\Delta`$. Write $`E=1.90(U+L)/\log L`$ and $`E_s=0.073(U+L)/(\log L)^2`$. At these large primes the square of the first component is $`(1+\chi_{547}\chi_{727})/2`$, and the same identity holds for B. The cross product is the mean of the four mixed pairs. The real two-by-two matrix $`G=\sum_{L\lt p\le U}J(p)e(p)e(p)^{\mathsf T}\log p`$ therefore has diagonal entries at least $`\Delta/8-E`$ and off-diagonal magnitude at most $`E`$. Here $`e=(e_A,e_B)`$. For each unit $`v`$, using $`2|v_Av_B|\le1`$ gives $`v^{\mathsf T}Gv\ge\Delta/8-2E\gt\Delta/20`$. Also $`|\sum J(v\cdot e)\log p|\le\sqrt2 E_s\lt\Delta/1000`$ and $`\sum J\log p\le\Delta/4+E_s\lt\Delta/3`$. All sums here are over primes in that interval.

Since $`|v\cdot e|\le\sqrt2\lt3/2`$, the positive and negative projected masses are separately greater than $`((\Delta/20)/(3/2)-\Delta/1000)/2\gt\Delta/64`$. Values of the chosen sign with magnitude at most $`1/32`$ contribute at most $`\Delta/96`$. Dividing the rest by $`3/2`$ gives log-prime mass greater than $`\Delta/288`$ for that sign above the stated cutoff. Taking $`U=4L/3`$ proves the prime bound. [Exact inequalities and all sign patterns](trace_moments.py)

The same matrix proves the full rank on the infinite range. Its nonzero prime vectors have four possible unoriented directions, the two axes and the two diagonals. If the directions represented by at least two distinct primes did not span the plane, take a unit vector perpendicular to them. The remaining singleton primes give quadratic mass at most $`8\log U`$. This contradicts the matrix bound whenever $`\Delta/20\gt8\log U`$. Thus there are two independent directions, each represented at least twice. In any common knot interval the two combined hats of an alleged relation are affine in $`p/y`$. For each repeated direction the corresponding affine projection vanishes at two different primes, hence throughout that interval. The independent directions force both hats to vanish.

For the three scales the distinct knots go from 1 to 8 with successive ratio $`\sqrt2`$. In a gap $`(c,d)`$ use $`L=(15c+d)y/16`$, $`U=(c+15d)y/16`$. These lie strictly inside the gap after scaling. Here $`(U+L)/\Delta=(8/7)(\sqrt2+1)/(\sqrt2-1)\lt7`$, $`\Delta\gt y/3`$ and $`L\ge y`$. For $`y\ge e^{366}`$ we have $`y/60\gt8\log(8y)`$ and the squares of all these primes exceed $`8y`$. Baker’s theorem and the repeated-direction argument now give the rank on every such algebraic start. More generally the same proof applies to a mixed list whenever its joint knot gaps contain subintervals meeting the two matrix inequalities, with squares above the largest supported norm. The prime chains in the exact transition proof connect the two ranges.

## The exact transition at the three scales

**Theorem.** Let $`y\gt0`$ be real algebraic. The six width-two scores from the two fields at $`y,\sqrt2 y,2y`$, each divided by $`\log2`$, are linearly independent over $`K`$ together with 1 **if and only if**

```math
y\gt\frac{79\sqrt2}{8}=13.9653589284\ldots .
```

For $`y\le79\sqrt2/8`$, the first two readings of B have upper norm endpoints at most 79. The only good prime norms in B through 79 are 49 from 7 and 79 from 79. At the upper endpoint the tent is zero, so 79 contributes to neither reading in this range. Both are algebraic multiples of $`\log7/\log2`$. The rank is deficient at the exact endpoint as well.

On the initial interval above that endpoint enumerate all rational primes through 12,000 and their prime-ideal norms on the two hosts. For each unramified prime, the common residue degree $`f`$ is the least of 1, 2 and 4 for which $`p^f\equiv\pm1\bmod16`$ and $`(\chi_u(p)\chi_v(p))^f=1`$ in that field. The $`8/f`$ ideals of norm $`p^f`$ have common sign $`\chi_u(p)^f`$. Their total logarithmic coefficient before the weight is $`-8\chi_u(p)^f\log p`$. Remove the common scalar, and remove each column’s nonzero algebraic square-root normalization. On each successive knot interval the entries of the resulting six-column coefficient matrix are affine functions of $`y`$, with coefficients in $`\mathbf Q(\sqrt2)`$.

The certificate selects six prime rows on each of **901** consecutive intervals. In the linear coordinate $`t`$ from its left to right endpoint the determinant has degree at most six. Its seven Bernstein coefficients have the same weak sign and the last coefficient has strict sign. The Bernstein basis is positive in the open interval, and its final term equals one at the right endpoint. Thus the determinant stays nonzero on the entire interval including its right endpoint. The replay reconstructs every knot in $`\mathbf Q(\sqrt2)`$, evaluates seven determinants by fraction-free elimination and inverts the seven-point interpolation matrix over the rationals. A single set of rational prime rows certifies an interval of real starts at once. Baker’s theorem gives the numerical rank at each algebraic one. This reaches 1,500 inclusive. [Exact minors](small_heights.py) · [Selected prime rows](evidence/initial-minors.tsv)

The six shorter prime chains above reach $`10^6`$. To continue from there to the analytic bound, take two increasing sequences of rational primes. A prime from the first has degree one only on A, and a prime from the second has degree one only on B. In each sequence the first prime is below one million, the last exceeds $`2^{531}`$, and each consecutive ratio is strictly less than $`\sqrt2`$. Thus each of $`(y,\sqrt2 y)`$, $`(\sqrt2 y,2y)`$ and $`(2y,2\sqrt2 y)`$ contains one prime of each type for every $`10^6\le y\le e^{366}`$. Indeed $`366\lt529\log2`$, giving $`2\sqrt2 e^{366}\lt2^{531}`$. Every chosen square is above $`8y`$. The diagonal hat in each band is positive and those for later scales vanish. These six rows again give two triangular blocks.

The table records 5,140 primes by exact integer offsets in five bins of every octave from 18 through 531. The 160 smaller primes are trial-divided. For the other 4,980 primes the table supplies $`p=1+k2^n`$ with odd $`0\lt k\lt2^n`$ and a witness $`a^{(p-1)/2}\equiv-1\bmod p`$. Every prime divisor of this odd $`p`$ would then have the order of $`a`$ divisible by $`2^n`$, making that divisor at least $`2^n+1\gt\sqrt p`$. This proves primality. The checker also verifies every gap, signed splitting pattern and finite omission. The analytic proof applies from $`e^{366}`$ onwards. Together the four overlapping ranges give the exact transition at every positive algebraic height. [Prime certificates](evidence/prime-bridges.tsv) · [Prime and gap verification](prime_bridges.py) · [Proth prime chains in prior computation](https://math.dartmouth.edu/~carlp/PDF/paper91.pdf)

## Exact crossings and deadlines

Let $`v_A,v_B`$ be real algebraic coefficients, at least one nonzero, and $`D_v=v_A D_A+v_B D_B`$. The existing unconditional theorem places $`[-5.85,5.85]^2`$ in the interior of the convex hull of actual moves after every $`H\ge1`$ and before one finite completion factor $`C`$. The supremum and the negative infimum of this projection both strictly exceed $`5.85(|v_A|+|v_B|)`$. Thus actual starts in the same eligible interval have $`D_v`$ on either side of any real algebraic $`t`$ with $`0\lt|t|\le5.85(|v_A|+|v_B|)`$. The projection is continuous and the interval is connected, so a single start hits $`t`$ exactly. Every such start is transcendental by the odd-prime-log argument. The same $`C`$ works for all these projections, all targets and both signs. For A alone the unconditional $`400H`$ cover through $`H=e^{12500}`$ gives both exact hits $`\pm1.9735`$ by that deadline. [Joint fixed-window theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/PROOF.md#one-fixed-window-for-all-directions) · [Explicit finite cover](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-waiting-bound-v1.3.0)

The two strict inequalities persist at nearby rational starts inside the eligible interval. The unique convex weight between those two projections that gives $`t`$ is transcendental. If it were algebraic, the nonzero algebraic $`t`$ would lie in the $`K`$-span of odd prime logarithms divided by $`a`$. The same obstruction applies to any finite number of algebraic starts with algebraic weights.

On each compact positive height interval only finitely many norms occur. Between consecutive breakpoints of the tents, a fixed algebraic combination of readings has the form $`\alpha\sqrt y+\beta/\sqrt y`$. The coefficients lie in the field generated by $`K`$ and the finitely many relevant prime logarithms, including $`a`$ for normalized scores. The equation for a threshold $`t\ne0`$ is $`\alpha w^2-tw+\beta=0`$ in $`w=\sqrt y`$. Its middle coefficient is nonzero, so there are finitely many solutions on each compact interval. For a nonzero algebraic target each solution $`y`$ is transcendental over $`\mathbf Q`$ and has degree at most two over this finite-log field, since $`w`$ does. Thus exact hits are locally finite despite occurring in every guaranteed multiplicative window.

Now fix algebraic $`C\gt2`$ and write $`k=C/2`$. For a nonzero real algebraic threshold and either sign, a component $`(l,u)`$ of the open qualifying set $`\{y:\sigma D_v(y)\gt t\}`$ serves exactly the starts $`H\in(l/k,u)`$. Indeed it must meet $`(H,kH]`$. The formula extends to components meeting zero or infinity. Near any finite positive $`H`$, the compact interval of eligible heights enlarged slightly past both ends meets only finitely many components, by local finiteness of the threshold roots. Every finite interior boundary of success or failure therefore comes from a finite positive $`u`$ or $`l/k`$ at such a component. Each is transcendental, since scaling a transcendental number by a nonzero algebraic number preserves transcendence. The boundary of a requirement for both directions is contained in the union of the individual boundaries. The first missed gain and first missed loss at factor 400 consequently have transcendental starting heights whenever they exist above the initial interval of success. This conclusion concerns $`H`$ itself. The numerical bounds in the directional result concern $`\log(2\sqrt2 H)`$ and assume its zero conditions, which GRH implies. [Directional first misses](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.3.0)

### How many rational starts does a vector target need?

In the guaranteed window for the two moves, every point of the closed square can be attained using at most **three rational starts** with real nonnegative weights. In the four-coordinate window both levels and both moves can be matched using at most **five rational starts**. Almost every target in the respective bodies requires that many, even when all positive algebraic starts are available. Thus the costs for rational starts are sharp. The antecedent common-window theorem supplies bounds of two and four starts when arbitrary real starts are allowed. [Common-window targets](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/PROOF.md#one-fixed-window-for-all-directions)

Indeed the image at rational starts is dense in the full continuous curve on the eligible interval. The closed convex hull of this rational image has every promised target in its interior. Such a target lies in the ordinary convex hull of the rational image. To see this, surround it by a small simplex in the interior, then approximate the simplex vertices by finite convex combinations of rational-start observations. The resulting simplex still contains the target in its interior, giving a finite representation. Compactness of the target body gives one finite pool of rational starts per window whose convex hull contains the entire body, by using a finite subcover of the interiors just constructed. Carathéodory’s theorem reduces each representation to $`d+1`$ of these observations in dimension $`d`$. For the lower bound, every tuple of at most $`d`$ algebraic starts spans a simplex of dimension at most $`d-1`$. There are only countably many tuples. Their union has $`d`$-dimensional Lebesgue measure zero. This lower bound concerns almost every target and does not select an individual algebraic target.

## Width refinement and exceptional real starts

For integer $`m\ge1`$ and $`b\gt1`$, define positive coefficients by

```math
\sum_{j=0}^{2m-2} c_jX^j
=(1+X+\cdots+X^{m-1})(b^{m-1}+b^{m-2}X+\cdots+X^{m-1}),
```

where each factor is 1 for $`m=1`$. At every real $`y,u\gt0`$ the identity is

```math
S_{i,b^m}(uy)=\frac{B_b}{B_{b^m}}\sum_{j=0}^{2m-2}c_j b^{j/2}S_{i,b}(b^juy).
```

On the geometric grid the narrow hat has jump polynomial $`(X-1)(X-b)`$ and the wide hat has jump polynomial $`(X^m-1)(X^m-b^m)`$. Their quotient is the displayed polynomial. The factors $`b^{j/2}`$ restore the normalization of the smaller tents. When $`b,u`$ are algebraic, the narrow columns are independent at every sufficiently large algebraic start. Thus the coefficients are unique over $`K`$ in this basis. [Knot insertion](https://www.sciencedirect.com/science/article/abs/pii/0010448580901542) · [Two exact calculations](verify_rank.py)

There is also a restriction on exceptional real starts for the rank theorem. Put $`\mathcal L=K(\log p:p\text{ rational prime})`$, a countable field. Fix a list of algebraic widths and scales and work strictly beyond its onset. For any coefficient vector over $`K`$ outside the asserted relation space, including a possible constant coefficient, multiply its value by $`\sqrt y`$. On each open prime-knot cell this gives a polynomial in $`\sqrt y`$ of degree at most two with coefficients in $`\mathcal L`$. It cannot vanish identically on a cell past the onset, since such a cell contains algebraic points where the rank formula holds. The breakpoints themselves are algebraic. Thus any exceptional real $`y`$ is transcendental over $`\mathbf Q`$ and has degree at most two over $`\mathcal L`$. There are at most countably many such starts. This describes where extra relations could occur and leaves their existence open.

## What an isolated missed deadline would force

Put $`L_{\rm odd}=K(\log p:p\text{ an odd prime})`$. Consider a nonzero real algebraic target $`t`$ for one signed move or for a fixed algebraic projection of the two moves. Its sign can be absorbed into the projection. Let $`C\gt2`$ be algebraic and $`k=C/2`$. Here a strict rational witness means a rational $`H`$ at which the entire closed interval $`[H,kH]`$ stays strictly below $`t`$.

**Theorem.** If a deadline fails at some real $`H\ge1`$ but has no strict rational witness with $`H\ge1`$, then $`(\log2)^2\in L_{\rm odd}`$. More precisely, a contact inside a missed window gives a degree-two homogeneous relation, or an isolated missed window gives a degree-four homogeneous relation, between finitely many prime logarithms. Both relations have a nonzero coefficient multiplying $`(\log2)^2`$. This implication is unconditional. Whether $`(\log2)^2`$ belongs to this odd-log field is open.

To prove it, write the signed projection with the common denominator $`a=\log2`$ cleared. Near a threshold hit its form is $`\alpha w+\beta/w`$ where $`w=\sqrt y`$, and $`\alpha,\beta`$ are $`K`$-linear forms in finitely many odd prime logarithms. Every threshold hit is transcendental, so it lies in the interior of a prime-knot cell. If it is a local maximum at the threshold, differentiation gives $`\alpha w^2=\beta`$ and hence

```math
t^2a^2=4\alpha\beta.
```

Both forms on the right are nonzero. If an algebraic $`H`$ misses the deadline, neither $`H`$ nor $`kH`$ hits the threshold. The value at $`H`$ cannot exceed the threshold either, since by continuity a strictly later start would then qualify. Thus either the closed interval stays strictly below the threshold, or equality is attained at an interior maximum. Every nontrivial interval of missed starts contains a rational start, so this covers that case.

Otherwise any missed start is locally isolated, by the local finiteness of component boundaries. If the only miss in the domain were $`H=1`$, the preceding rational argument would already apply. At an isolated missed $`H\gt1`$ a qualifying component must end at $`H`$ and another must start at $`kH`$. Indeed the open success intervals $`(l/k,u)`$ cover arbitrarily close starts on the left and right, and locally only finitely many are relevant. At both endpoints the signed move equals $`t`$. Immediately to the left of $`H`$ the move qualifies, while immediately to the left of $`kH`$ it cannot qualify.

On neighbourhoods of $`y=H`$ write the two cleared moves at $`y`$ and $`ky`$ as $`\alpha_1w+\beta_1/w`$ and $`\alpha_2w+\beta_2/w`$. Scaling by algebraic $`k`$ preserves the odd-log coefficient property. The two pairs of coefficients cannot coincide because the qualifying signs on the left differ. Subtracting the two threshold equations gives $`Pw^2=Q`$, with $`P=\alpha_1-\alpha_2`$, $`Q=\beta_2-\beta_1`$. Both are nonzero. Substitution then gives the explicit identity

```math
t^2a^2(\alpha_1-\alpha_2)(\beta_2-\beta_1)=(\alpha_1\beta_2-\alpha_2\beta_1)^2.
```

This proves the theorem, including singleton failures where equality at a boundary would escape a rational starting-height search. For a joint requirement for both signs, apply the same argument to whichever sign fails. [Quadratic elimination in Lean](AlgebraicHeights.lean)

Under the additional condition $`a^2\notin L_{\rm odd}`$, every actual failure at an algebraic deadline therefore has a strict rational witness somewhere. A dovetailed search can enumerate rational $`H\ge1`$ and positive precisions, enumerate the finitely many norms in the interval and enclose the projection on closed subdivisions. Continuity and the positive margin guarantee termination at a genuine strict witness if a failure exists. This is a completeness criterion for a search, with no practical runtime bound. The unconditional finite alternative above is available without the added independence condition. Schanuel’s conjecture would imply that condition. Indeed any finite list of distinct prime logarithms is rationally independent and its exponentials are rational, so under the conjecture the logarithms are algebraically independent. [Statement and context of Schanuel’s conjecture](https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/ColloquiumDeGiorgiSchanuelConjecture2014.pdf)
