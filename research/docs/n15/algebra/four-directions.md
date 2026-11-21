# 用四个方向消去 N=15 覆盖问题的接触链

## 已验证的半径

2026-11-18 完成重放。

$`r_D(15)`$ 的严格十进制界限

- 下界为 0.318142930859262836359491181649935005186663866157651845538266。
- 上界为 0.318142930859262836359491181649935005186663866157651845538267。

精确对象为既有 148 元整系数系统在指定有理盒中的唯一根 $`\xi`$。从零编号，$`r_D(15)=\sqrt{\xi_{76}}`$，第 $`i+1`$ 个圆心是 $`(\xi_{20+2i},\xi_{21+2i})`$，其中 $`0\le i<15`$。系统、根盒中心与预条件器实际列在 [system15.py](../../../n15/proofs/global-final/replay/core/system15.py) 和 [root15.json](../../../n15/proofs/global-final/replay/core/root15.json)。盒中心为文件中的 xnum/Qx，逐坐标半宽可取 $`10^{-80}`$。

严格压缩验证给出盒内唯一性。连续覆盖证明使用 39 个三角形与 11 个圆弓。任意更小半径覆盖经一般位置扰动、极小子覆盖和 Voronoi 对偶归约为五族辅助杆图。独立枚举审计覆盖 11,950,884 个类型，度量不等式排除 11,907,870 个。剩余 43,013 个非候选类型和一个候选类型的分支逐叶检查完毕。本次在空的临时目录用固定版源代码重新编译、验算，五份归档代表表均重新经过独立审计。另在独立目录重放了 $`148\to56\to29\to27\to25`$ 元根的代数连接，逐字节核对最终的轴根与本次消元输入。人读论证、来源和两份实际报告分别见[全局证明](../proofs/global-final/global-proof.zh.md)、[固定版上游证明包](https://github.com/VonEquinox/DiskCoveringSolve/tree/089584973f41e6374ce110f6e5b5ee456df19352/cover15/proof_bundle)、[全局报告](../../../n15/proofs/global-final/master-20261118.json)和[根链报告](../../../n15/proofs/global-final/axis-chain-20261118.json)。

上面的十进制区间中恰有最优半径的 2,869 次最小多项式的一个实根。[完整解答](../solution.zh.md)给出全部整系数和有序圆心。[半径多项式核验](radius-identity.md)证明该根与原构型相合。下节另给一种利用原 148 元系统直接隔离最优半径所有共轭根的方法。

## 用有理端点隔离半径的所有共轭根

记 $`f`$ 为上述 148 元系统，$`q`$ 为 `root15.json` 中的精确盒中心。令 $`z_0=q`$，$`u_0`$ 为上述严格十进制下界。下式的 768 次复合用有限次四则运算指定两个有理数。矩阵逆也可按伴随矩阵除以行列式理解。

```math
\begin{aligned}
z_{j+1}&=z_j-f'(z_j)^{-1}f(z_j),\\
u_{j+1}&=\frac{u_j^2+(z_{j+1})_{76}}{2u_j},\\
a&=u_{768}-2^{-10^{223}},\\
b&=u_{768}+2^{-10^{223}}.
\end{aligned}
```

**根分离。** 给原系统加上 $`u^2-z_{76}=0`$，得到 $`n=149`$ 条整系数二次方程 $`f_i=0`$。原根非奇异且 $`u=\alpha>0`$，故增广根也非奇异。各式系数绝对值之和至多 1,024。取独立参数 $`v`$，考察 $`x_i^3=v f_i(x)`$。对任一次数相容单项式序，其首项两两互素并且首一，标准单项式为每个指数小于 3 的 $`d=3^{149}`$ 个单项式。也可逐一替换立方；不同坐标的替换在共同的乘积上可归并。用这组基表示乘以最后一个坐标 $`u`$ 的矩阵。各列起始单项式的总次至多 299，每次替换总次至少下降一。每个矩阵元素作为 $`v`$ 的整系数多项式，其系数绝对值之和至多 $`1024^{299}`$。

该矩阵特征多项式 $`H(v,t)`$ 关于 $`t`$ 首一，其系数绝对值之和至多 $`d!(1+1024^{299})^d<2^{10^{76}}`$。令 $`L(t)`$ 为 $`v`$ 的最高非零次项的系数。把变形写为 $`f_i(x)=\epsilon x_i^3`$，非奇异性给出 $`\epsilon=0`$ 附近通过目标根的解析解。在 $`\epsilon\ne0`$ 时，该解的最后一个坐标使 $`H(1/\epsilon,t)`$ 为零。乘以 $`\epsilon^{\deg_v H}`$ 并取极限得 $`L(\alpha)=0`$。因而 $`L`$ 是非零整系数多项式，次数至多 $`d`$，高度至多 $`2^{10^{76}}`$。

设 $`P`$ 是 $`\alpha`$ 的本原最小多项式，则 $`P`$ 在整数环整除 $`L`$，其首项绝对值至多 $`H_0=2^{10^{76}}`$。它的每个根按 $`L`$ 的 Cauchy 界至多有模 $`1+H_0`$。若 $`P`$ 有两个不同根，取其间距 $`s`$，在非零整数判别式的根差乘积中单独保留这对。由

```math
\begin{aligned}
1\le|\mathrm{Disc}P|\\
&\le H_0^{2d-2}s^2\\
&\quad{}\cdot[2(1+H_0)]^{d(d-1)-2}
\end{aligned}
```

和 $`d<10^{72}`$ 得 $`s>2^{-10^{222}}`$。一次多项式自动满足所需隔离结论。

**端点有效性。** 原根箱验证在半宽 $`10^{-60}`$ 的盒内给出 $`\|I-Yf'(z)\|_\infty<1/2`$，固定预条件器满足 $`\|Y\|_\infty<152`$。因此盒内每个 $`f'(z)`$ 可逆且逆范数小于 304。二次方程的 Taylor 余项至多 $`1024\|z-\xi\|_\infty^2`$。原 148 坐标的初始误差小于 $`10^{-80}`$，故之后的 Newton 点一直在大盒内。平方根坐标有

```math
u_{j+1}-\alpha=\frac{(u_j-\alpha)^2+(z_{j+1})_{76}-\xi_{76}}{2u_j}.
```

取两部分误差的最大值 $`e_j`$，初值小于 $`10^{-60}`$。只要 $`u_j>1/4`$ 就有 $`e_{j+1}<10^{10} e_j^2`$，此界也归纳保证 $`u_j>1/4`$。迭代给出 $`e_{768}<2^{-2^{768}}<2^{-10^{223}}`$。所以 $`0<a<\alpha<b`$，区间全长小于上面的根间距，$`P(a)P(b)\ne0`$，且其中只有实根 $`\alpha`$。端点以有理算式交付，未展开其分子与分母。实际证书检查方程系数、两个源根箱、预条件器和所有整数指数比较。[常数核验程序](../../../n15/algebra/four-directions/conjugate_interval.py) · [实际核验记录](../../../n15/algebra/four-directions/conjugate-interval.json)

## 四方向参数化

沿用[八变量系统](contact-transfer-20261024/proof.zh.md)的 $`u,v,w,P=(p_x,p_y),z,h,y,x_6,r`$。用 $`A,B,C,D`$ 表示以下四个单位向量，后文也把它们视为单位复数。此处 $`D`$ 表示方向变量。

```math
\begin{aligned}
rA&=P-(u,0),\\
rB&=(v,w)-P,\\
rC&=((z-p_x)/2,h-p_y),\\
rD&=(x_6-1,y).
\end{aligned}
```

每个方向的独立实参数是 $`s=Y/(1+X)`$，逆为 $`(X,Y)=((1-s^2)/(1+s^2),2s/(1+s^2))`$。设 $`k=A_x+2C_x+2(D_x-B_x)-1`$。由 $`u=1+r(2D_x-2B_x-1)`$ 及 $`(z,p_y)=(1+rk,rA_y)`$ 在单位圆上得到

```math
r=\frac{-2k}{k^2+A_y^2}.
```

以此顺序求 $`u,P,(v,w),z,h,y`$，并令 $`t=2((v-u)p_y-w(p_x-u))`$，即得八变量系统的四条载体方程。反过来，这四条方程使上述四个向量的模均为 1。例如首式 $`t^2=L(4r^2-L)`$，其中 $`L=(v-u)^2+w^2`$，给出 $`|P-(u,0)|^2=|P-(v,w)|^2=r^2`$。其余三式给出 $`z^2+p_y^2=1`$ 及方向 $`C,D`$ 的模。恢复公式和全部除法在源根盒上均可执行。这给出目标邻域的双向有理变换。

边界接点 $`q`$ 的下一圆心若为 $`q+rV`$，另一边界接点为

```math
M_V(q)=(q+rV)/(1+rV^{-1}q).
```

这是将 $`q`$ 关于过原点和圆心的直线反射。左链从 $`E=z+i p_y`$ 起步。三个连续圆心分别是 $`E-r/C,\ q_1+rB,\ q_2-rC`$。第一式代入定义即得，后两式来自两组等半径圆交点构成平行四边形的向量恒等式。因此左端接点为 $`\ell=M_{-C}M_B M_{-1/C}(E)`$，右端接点为 $`V=M_D(1)`$。另一个内接点为

```math
W=1+r(D+D^{-1}-1+A).
```

具体地，旧系统右端所选圆的圆心是 $`(p_x+2rB_x,p_y)`$。在旧右端圆方程中代入 $`v-u=r(A_x+B_x),\,w=r(A_y+B_y)`$ 及两条单位模恒等式即得此坐标。目标的第 5 号圆心（从零编号）到 $`\ell,V,W`$ 等距 $`r`$。另外两个边界接点到该圆心的平方距离都严格大于 $`r^2`$，且 $`\ell,V,W`$ 的有向面积非零。这些严格不等式在邻域保持，三点外心唯一恢复该圆心。因而所选分支的闭环可写为三角形的外接半径等于 $`r`$。四变量的旧方程拉回与这条新闭环的临界方程均通过独立的有理盒压缩。见[坐标及几何公式](../../../n15/algebra/four-directions/chart.py)、[四变量验算](../../../n15/algebra/four-directions/verify_chart.py)和[根证书](../../../n15/algebra/four-directions/target.json)。

## 消去中间两杆

用复代数坐标时，单位方向的共轭写成逆元，方程的系数属于 $`\mathbb Q`$。现在以 $`E,A,D`$ 为三个自由方向，设

```math
\begin{aligned}
r&=\frac{E-E^{-1}}{A-A^{-1}},\\
T&=\frac{E-1}{r}-A-D-D^{-1}+1\\
 &=C+C^{-1}-B-B^{-1},\\
U&=1-r^2.
\end{aligned}
```

固定 $`E,A`$ 即固定半径。三个分式线性映射的矩阵乘积为

```math
\begin{pmatrix}U+r^2CT&-r(T+U/B)\\-r(T+UB)&U+r^2T/C\end{pmatrix}.
```

要求此矩阵把 $`E`$ 送到 $`\ell`$ 给出 $`aC+b/C+cB+d/B+k=0`$，其中

```math
\begin{aligned}
a&=-r^2TE,\\
b&=\ell r^2T,\\
c&=-\ell rUE,\\
d&=rU,\\
k&=\ell(U-rTE)-UE+rT.
\end{aligned}
```

合并 $`C+C^{-1}=B+B^{-1}+T`$。记 $`x=b+c,\ y=b+d,\ z_0=a-b,\ w_0=k+bT`$，并令 $`p(X)=xX^2+w_0X+y`$。消去 $`C`$ 后得到

```math
\begin{aligned}
Q(X)&=p(X)\big[p(X)+z_0(X^2+TX+1)\big]\\
&\quad{}+z_0^2X^2,\\
Q(B)&=0.
\end{aligned}
```

这是从两条线性等式解出 $`C,C^{-1}`$，再将它们相乘等于 1 的结果。目标处 $`z_0\ne0`$。令 $`K=W-V,\ J=V\overline W-1`$。用 $`\bar\ell=1/\ell,\bar V=1/V`$ 代入三点外接半径公式并约去非零弦长，得到二次式

```math
H(\ell)=(W-\ell)(\ell\overline W-1)KJ-r^2(K+\ell J)^2.
```

在 $`E,A,D`$ 固定时可局部以 $`B`$ 为参数，因为 $`C-C^{-1}\ne0`$。固定三者的变化保持 $`r`$ 不变。目标闭环及 $`r`$ 的微分成比例，故沿该方向 $`H_\ell\ell_B=0`$。证书给出 $`H_\ell\ne0`$，进而 $`\ell_B=0,\ Q_B=0`$。另外 $`Q_{BB},Q_\ell`$ 和四次首项均非零。

将四次判别式中的 $`z_0^4`$ 剥去，记所得多项式为 $`\Delta(\ell)`$。它可用一个三次判别式计算。其三次多项式按升幂的系数为

```math
\begin{aligned}
&z_0^2,\\
&Tw_0z_0+w_0^2-4xy-2(x+y)z_0+z_0^2,\\
&T^2xy-Tw_0(x+y)+w_0^2+(x-y)^2\\
&\qquad{}-4xy-2(x+y)z_0,\\
&(x-y)^2.
\end{aligned}
```

对独立符号 $`x,y,z_0,w_0,T`$ 的整系数展开验证了 $`\mathrm{Disc}_X Q=z_0^4\Delta`$。因此可在 $`H`$ 的二次商代数取范数以消去 $`\ell`$。实际盒内 $`H`$ 的另一根上原四次判别式非零，$`H(-E)\ne0`$。选定接点处判别式对 $`\ell`$ 的导数也非零。见[消元公式](../../../n15/algebra/four-directions/eliminate.py)、[验算程序](../../../n15/algebra/four-directions/verify_elimination.py)和[三变量证书](../../../n15/algebra/four-directions/projection.json)。

## 一个 104 次的方向多项式

令 $`R=\mathrm{Res}_{\ell}(H,\Delta)`$，采用次数 2、8 的齐次结式，因而 $`R=h_2^8\mathrm{Norm}_{H/h_2}\Delta`$。设

```math
F(D)=D^{36}(D+r)^{16}R.
```

**次数界。** 固定 $`E,A,r`$，其余系数属于特征零有理函数域。只有 $`D=0,-r`$ 可能是 $`R`$ 的有限极点。$`H`$ 的各项在 $`-r`$ 的极点阶数至多 2，贡献至多 16。对于另外两处，设 $`v=1/T`$，把 $`\Delta`$ 写成接点齐次坐标 $`(j:k)`$ 的八次式。程序逐项展开 $`v^{22}\Delta`$，其最低 $`v`$ 次数在通用 $`(j,k)`$ 为 10，在 $`(r,-v+v^2s)`$、$`(-v+v^2s,r)`$ 均为 14。这三个等式在独立参数上成立。

为看清 $`H`$ 的两个根，设 $`d=D,\ U_0=1+r(A-1),\ U_1=1+r(A^{-1}-1)`$。直接相减得

```math
\begin{aligned}
K&=\frac{(U_0+r^2-1)d^2+r(1+U_0)d+r^2}{d(d+r)},\\
J&=\frac{d}{d+r}\bigl[r^2d^2+r(1+U_1)d\\
&\quad{}+U_1+r^2-1\bigr].
\end{aligned}
```

当 $`d\to0`$，$`d^2H`$ 正则，极限是一条有限根与一条简单无穷根。后者满足 $`1/\ell=d/r+O(d^2)=-v/r+O(v^2)`$。当 $`d\to\infty`$，$`d^{-4}H`$ 正则，有一条简单根 $`\ell=1/(rd)+O(d^{-2})=-v/r+O(v^2)`$。这些由上式代入 $`H`$ 的常数项、一次项及二次项相除即得，均在通用系数域计算。两处各有一根使齐次判别式的极点阶数由 12 降到至多 8。由结式的齐次乘积公式，$`R`$ 在零处极点阶数至多 $`2\cdot8+12+8=36`$，在无穷处增长次数至多 $`4\cdot8+12+8=52`$。于是 $`F`$ 是次数至多 $`36+16+52=104`$ 的多项式。通用恒等式在清分母后可特化至目标。

在 $`E=2,A=3`$，特征零展开实际得到 104 次。去内容后与导数的整数结式有 17,866 位十进制数字，模 65,521 等于 59,522，证明该截面无重根。全部截面系数、结式摘要与目标处非零因子见[二参数证书](../../../n15/algebra/four-directions/binary.json)和[验算程序](../../../n15/algebra/four-directions/binary_projection.py)。

## 实际半径在平面曲线上的分支

在实半角参数中，取未剥去四次幂的二次范数之实部为 $`k`$。在实际盒内认证 $`k_{s_E}\ne0,\ k_{s_Ds_D}\ne0`$。用 $`Q_\ell,Q_{BB},H_\ell\ne0`$ 可知原几何闭环沿自由 $`B`$ 方向为非退化二阶接触。先用隐函数定理解闭环对 $`B`$ 的驻定方程，再代回闭环，得到一张光滑实曲面。其梯度非零来自原闭环光滑且沿 $`B`$ 的导数为零。曲面上 $`k=0`$，且 $`k`$ 的梯度非零，故在目标邻域两张曲面重合。被约去的因子在此邻域均非零。

在曲面上 $`r`$ 仍在目标处驻定。$`r`$ 仅依赖 $`E,A`$，且局部对 $`s_E`$ 的偏导非零。由此得 $`k_{s_D}=0`$。三式 $`k,\ k_{s_D},\ k_{s_E}r_{s_A}-k_{s_A}r_{s_E}`$ 的独立压缩通过。沿 $`k=k_{s_D}=0`$，可用 $`s_A`$ 解出 $`s_E,s_D`$，因为对应雅可比行列式为 $`k_{s_E}k_{s_Ds_D}\ne0`$。在此实曲线上，$`D`$ 方向与曲面相切，故复数方程也有 $`F=F_D=0`$。沿曲线半径在目标处的一阶导数为零。若它局部恒定，则会给出上述三式的非孤立零点，与压缩矛盾。

写 $`x=r`$，只保留方向 $`A`$。辅助方向 $`E`$ 满足 $`E^2-x(A-A^{-1})E-1=0`$，另一根为 $`-1/E`$。将 $`F`$ 中的半径保留为独立变量 $`x`$，令

```math
\begin{aligned}
Z(x,A)&=A^{14904}\mathrm{Norm}_E\bigl[\\
&\qquad{}\mathrm{Res}_D(x^{24}F,\partial_D(x^{24}F))\bigr].
\end{aligned}
```

这是整系数平面多项式。可先用 $`A^2D^2(D+x)^2H`$ 和 $`(xD)^{12}\Delta`$ 计算第一个结式，再约去 $`D^4`$ 及坐标单项式。固定 $`x`$ 时，$`H`$ 中 $`A`$ 的指数在 $`[-2,2]`$。$`\Delta`$ 对 $`T`$ 最高 12 次，显含 $`E`$ 至多 8 次。对 $`F`$ 的系数，$`A^iE^j`$ 满足 $`-16\le i\le40,\ 0\le j\le40,\ i+j\le56`$。$`x^{24}F`$ 对 $`x`$ 至多 120 次，因为清分母后的 $`H,\Delta`$ 对 $`x`$ 分别至多 6、36 次。再取一次 104 与 103 次的结式，将这些指数界乘以 207。

当 $`A\to\infty`$，两根 $`E`$ 分别为 $`O(A)`$ 和 $`O(A^{-1})`$，范数增长次数至多 $`(56+40)\cdot207`$。当 $`A\to0`$，两根分别为 $`O(A^{-1})`$ 和 $`O(A)`$，总极点阶数至多 $`(56+16)\cdot207=14904`$。对 $`x`$，两个单项式的范数交叉项可用 $`EE'=-1`$ 和 $`E+E'=x(A-A^{-1})`$ 的幂和递推化简，附加次数至多 $`40\cdot207`$。故 $`\deg_A Z\le34776,\ \deg_x Z\le57960`$。

**目标因子选择。** 将四变量盒进一步收缩到半宽 $`10^{-350}`$，在 128 次单位根上用复球算术求值并反变换，包住实际 $`F`$ 的全部系数。两个 $`E`$ 分支的首项均避零。已证明目标方向 $`D_*`$ 为重根，可在目标分支作两次有误差外包的综合除法。余下的 102 个根全部分别隔离，且余因子在 $`D_*`$ 避零。另一分支的 104 个根也都分别隔离，两支没有公共根。见[实际纤维报告](../../../n15/algebra/four-directions/target-fiber.json)和[复验程序](../../../n15/algebra/four-directions/target_fiber.py)。

因此 $`F`$ 在目标恰有一次二根相合。因为 $`F_{DD}\ne0`$，该临界方向能由 $`F_D=0`$ 在邻域唯一延拓，而 $`F`$ 沿横截方向的一阶导数非零来自上述曲面的梯度。局部用首一二次因子分出这对根，与其余简单根的结式均为单位。该二次的判别式横截简单消失。这说明 $`Z`$ 在目标的局部因子重数为 1。变换 $`(s_E,s_A)\mapsto(x,A)`$ 在这里可逆，因为 $`A\ne\pm1`$、$`E+1/E\ne0`$。因而这项结论也适用于 $`(x,A)`$ 坐标。

在 $`\mathbb Q[x,A]`$ 对 $`Z`$ 作平方自由分解，令 $`S_1`$ 为重数恰为一且依赖 $`A`$ 的不可约因子之积。目标落在其中一条曲线上。前述沿 $`s_A`$ 的驻定性，及 $`dA/ds_A\ne0`$，给出 $`S_1=\partial_A S_1=0`$。故

```math
G(x)=\mathrm{Res}_A(S_1,\partial_A S_1)
```

是一个**非零且在实际半径消失**的单变量关系的代数表达式。非零来自 $`S_1`$ 在 $`\mathbb Q(x)[A]`$ 平方自由。以给出的通用次数界可得 $`\deg G\le4\,031\,175\,960`$。此界只作上界。后续计算直接重建并验证包含目标的不可约因子，避免展开整个 $`G`$。

## 全参数整系数输入

直接在五变量整数环消去 $`\ell`$ 的试算在进程占用接近 6 GB 的作业上限时退出，返回 Windows 代码 0xC0000005。随后改为先在 $`w=AE`$ 的二次代数中规约，关系为 $`w^2=r(A^2-1)w+A^2`$。对 $`H=h_2\ell^2+h_1\ell+h_0`$ 的两根，用 $`t_0=2,\ t_1=-h_1,\ t_k=-h_1t_{k-1}-h_0h_2t_{k-2}`$ 计算齐次幂和。若待取范数的多项式为 $`\sum_{i=0}^m p_i\ell^i`$，其齐次结式为

```math
\sum_i p_i^2 h_0^i h_2^{m-i}
+\sum_{i\lt j} p_i p_j h_0^i h_2^{m-j}t_{j-i}.
```

这是将 $`h_2^mP(\ell_1)P(\ell_2)`$ 的对角项与一对交叉项合并所得。[完整整系数表](../../../n15/algebra/four-directions/carrier.csv.gz)实际列出了

```math
P(r,A,D)+wQ(r,A,D)=A^{48}r^{24}F(D).
```

文件首行为 JSON，之后每行为分量编号、$`r,A,D`$ 的非负指数以及有符号整数系数。$`P,Q`$ 分别有 315,019 和 329,394 项。两者都含 $`r^{48}A^{12}`$，可在目标处约去。已逐项重建全部系数，并在三个有理方向组合上独立重算全部 105 个 $`D`$ 系数。[生成器](../../../n15/algebra/four-directions/structured_expansion.py) · [验算器](../../../n15/algebra/four-directions/verify_carrier.py) · [实际记录](../../../n15/algebra/four-directions/carrier-verified.json)

在 $`E=2,A=3`$，取 104 次方向多项式的本原整系数表示并模素数 1093。其次数保持，$`D^{1093^{104}}\equiv D`$，且它与 $`D^{1093^k}-D`$ 在 $`k=8,52`$ 时的最大公因子均为 1。由 $`104=2^3\cdot13`$ 及 Rabin 判据，该截面在有理数域不可约。清去参数分母后，首项在这个特化处非零，故通用方向多项式在 $`\mathbb Q(E,A)[D]`$ 也不可约。这个不可约性结论的变量是 $`D`$。

对 $`w`$ 取范数所得 $`N=P(P+r(A^2-1)Q)-A^2Q^2`$ 有 1,845,061 项，三变量次数依次为 248、148、208。在实际根处，另一支和目标支的所有根盒互不相交，故这个范数也保留目标处唯一的二根相合。取范数后的完整特征零 $`D`$ 结式在 300 秒受限试算内没有输出。这两次试算的停止原因分别见[直接展开记录](../../../n15/algebra/four-directions/expansion-run.json)和[结构化展开记录](../../../n15/algebra/four-directions/structured-run.json)。[特征零展开程序](../../../n15/algebra/four-directions/expand_radius.py)

## 实算与复现

在 $`\mathbb F_{998244353}`$ 上计算四个完整截面。固定 $`x`$ 时，从报告中恢复 $`Z`$ 须再乘非零常数 $`x^{9936}`$。每个截面按上述次数界插值，另用八点校验，并用独立的单变量有理函数实现直接复核六组结式，再逐因子相乘。截面的平方自由分解在特化后进行。

第二张表列出重数一层首一化后与其导数的模结式。

| 特化 | 变量 | 截面次数 | 单重层次数 |
| --- | ---: | ---: | ---: |
| $`x=9/16`$ | $`A`$ | 20,888 | 1,320 |
| $`x=2/5`$ | $`A`$ | 20,888 | 1,320 |
| $`A=3`$ | $`x`$ | 49,296 | 1,416 |
| $`A=1/3`$ | $`x`$ | 49,296 | 1,416 |

| 特化 | 模结式 |
| --- | ---: |
| $`x=9/16`$ | 907706728 |
| $`x=2/5`$ | 170316187 |
| $`A=3`$ | 136466225 |
| $`A=1/3`$ | 202346956 |

这些是确切的截面计算。特化可能改变因子及其重数，表中的数值不能确定最优半径的最小多项式次数。原始系数和验证参数保存在[第一个半径截面](../../../n15/algebra/four-directions/slice-9-16.json.gz)、[第二个半径截面](../../../n15/algebra/four-directions/slice-2-5.json.gz)、[方向截面](../../../n15/algebra/four-directions/slice-a3.json.gz)及[互反方向截面](../../../n15/algebra/four-directions/slice-a-third.json.gz)。[Python 驱动](../../../n15/algebra/four-directions/radius_slice.py)和[C++ 求值器](../../../n15/algebra/four-directions/radius_slice.cpp)同时交付。

### 利用方向互反缩小模计算

沿实际临界曲线，$`r`$ 为实数且 $`\bar A=1/A`$。因 $`S_1`$ 的系数为有理数，$`S_1(r,A)=S_1(r,1/A)=0`$ 在整条实解析弧上成立。因此清去单项式后，两式必须有包含该弧的公共代数分支。非固定的 $`A`$ 满足 $`A\ne\pm1`$，在这个邻域可用 $`X=A+A^{-1}`$ 作方向参数。对偶次首一互反多项式 $`f`$，用 $`A^k+A^{-k}=X(A^{k-1}+A^{1-k})-(A^{k-2}+A^{2-k})`$ 写出 $`A^{-\deg(f)/2}f(A)`$ 的 $`X`$ 多项式。

在上述两个固定半径上，1,320 次层与其倒序多项式的首一最大公因子均为 1,312 次，换元后均为 656 次。在第二种模算法中，先从整系数表同时约掉 $`r^{48}A^{12}`$，随后 $`P,Q`$ 对 $`A`$ 至多 66 次。取范数时 $`w=AE`$ 在 $`A\to\infty`$ 的两支分别为 $`O(A^2)`$ 和 $`O(1)`$。固定半径处的方向结式范数因而对 $`A`$ 至多 $`(68+66)\cdot207=27738`$ 次。每个样本在有 65,536 个单位根的网格上选择 27,739 个分裂点，逐点计算两支的结式并插值，另取四点检查。前两个有理半径的 1,321 个系数均与早先独立程序逐项相同。[系数表求值器](../../../n15/algebra/four-directions/plane_samples.cpp) · [样本驱动](../../../n15/algebra/four-directions/sample_plane.py)

### 从方向曲线到模半径关系

在整数半径 2 到 1503 中有 1,500 个截面的重数一层次数恰为 1,320。固定 $`A=3`$ 的半径截面用来统一每个方向多项式的标量，以 1,480 个半径插值，另外 20 个半径核对每个系数。所得 $`(r,A)`$ 多项式的实际双次数为 $`(1436,1320)`$。对互反公因子，改用 $`X=10/3`$，以 $`A=3,1/3`$ 两个截面的首一公因子统一标量。1,470 个半径给出双次数 $`(1424,656)`$ 的 $`C(r,X)`$。另外 30 个半径以及两个独立有理半径逐系数相符。截面的次数和全参数因子的次数可以不同。这是经留出点验证的有限域插值，尚无通用恒等式把整条 $`C`$ 与特征零目标曲线等同起来。[全部原始样本](../../../n15/algebra/four-directions/plane-samples.zip) · [插值程序](../../../n15/algebra/four-directions/fit_plane.py) · [对称曲线与存储格式](../../../n15/algebra/four-directions/plane-symmetric.json)

对于已列出的 $`C`$，末级结式 $`J(r)=\mathrm{Res}_X(C,C_X)`$ 的次数至多 $`1424(2\cdot656-1)=1,866,864`$。在全部 $`2^{21}`$ 个单位根上计算齐次结式并反变换，得实际次数 1,759,328。独立地在 $`r=9/16,2/5,3/7`$ 直接取结式，三值逐一相符。平方自由部分有 694,039 次，其中在 $`J`$ 中重数恰为一的因子之积有 3,152 次。这一步对给定模曲线精确。[末级程序](../../../n15/algebra/four-directions/project_mod_radius.py) · [次数、系数文件以及交叉检查](../../../n15/algebra/four-directions/mod-radius.json)

逐个将 3,152 次层的不可约因子 $`f`$ 代入 $`C`$，在 $`K=\mathbb F_p[r]/(f)`$ 中检查 $`\gcd(C,C_X)`$。有九个因子保持方向次数并给出有限普通二重根 $`X`$，次数是 1、2、12、21、64、253、655、776、1085，合计 2,869。其余因子处 $`C`$ 的方向首项消失。对这九组再以 $`A^2-XA+1=0`$ 对 208 次的方向方程取范数，方向重根的公因子均为 $`D^2+gD+1`$。在这个二次代数中，由范数的一次余式求出 $`A`$ 并用 $`E=-P/(AQ)`$ 恢复 $`E`$。直接验证原特征零公式约化模 $`p`$ 后的 $`F=F_D=F_E r_A-F_A r_E=0`$ 与雅可比可逆性，九个因子全部通过。这给出对应临界点自身的模证书。这一步取得了模坐标；其特征零连接由后续整数恒等式给出。[半径与第一方向](../../../n15/algebra/four-directions/mod-fibers.json) · [原临界方程回代](../../../n15/algebra/four-directions/mod-directions.json)

### 用原方程提升半径系数

在每个非奇异模纤维内，先把八次分支式约化模接点二次式，其一次余式解出左接点 $`\ell`$。在 $`\ell`$ 处，四次式 $`Q`$ 与 $`Q_B`$ 的一次公因子解出 $`B`$。保留五个变量 $`(E,A,D,\ell,B)`$。记 $`T_i=H_iQ_\ell-Q_iH_\ell`$，其中偏导固定其余四个变量。五条方程为

```math
H=Q=Q_B=T_D=T_Er_A-T_A r_E=0.
```

九个剩余域都通过这五式的原式回代。在模 $`p^k`$ 的一支，以 $`R`$ 表示半径，$`M(R)`$ 为首一半径多项式，有向坐标的基为 $`1,D_0`$，满足 $`D_0^2+g(R)D_0+1=0`$。把原方程模 $`p^{2k}`$ 求值，所有残差被 $`p^k`$ 整除。除去该因子后，用旧精度的雅可比求修正量。新半径为 $`R+\delta(R)`$，其中 $`\delta^2=0\pmod {p^{2k}}`$。在新坐标中 $`M`$ 改为 $`M-(M'\delta\bmod M)`$。对系数多项式用相同的一阶换元变换。再将基中的第二个元素换成新有向根，使其和逆元之和无有向分量。每次换元保持半径多项式首一。

用三元判别式方程单独提升到 $`p^2`$，再用五元接触方程从同一模纤维提升到 $`p^8`$，两种计算的全部 2,870 个半径系数在 $`p^2`$ 上逐项相等。已独立检查每支的方程与半径坐标。各精度的乘积给出九支全部 2,869 个模半径。从 $`p^{256}`$ 重建的本原整系数及其与实际半径的连接见[特征零核验](radius-identity.md)。[原始接点](../../../n15/algebra/four-directions/mod-anchors.json) · [代数提升程序](../../../n15/algebra/four-directions/padic_radius.py) · [三元提升](../../../n15/algebra/four-directions/padic-p2.json.gz) · [五元提升](../../../n15/algebra/four-directions/padic-p64.json.gz)

从仓库根目录运行以下命令。前两条命令使用固定 Git 归档，报告输出路径必须尚不存在。第三条命令重新计算八份代数报告并逐字段比较。

~~~powershell
& scripts/serialized_compute.ps1 python -X utf8 -B `
    n15/proofs/global-final/replay_archived_inputs.py --output global-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B `
    n15/proofs/global-final/replay_archived_inputs.py --axis-chain --output `
    axis-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B n15/algebra/four-directions/replay.py
& scripts/serialized_compute.ps1 python -X utf8 -B `
    n15/algebra/four-directions/verify_carrier.py --rebuild --output `
    carrier-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B `
    n15/algebra/four-directions/conjugate_interval.py --output interval-fresh.json
$dir = 'n15/algebra/four-directions'
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/fit_plane.py" `
    "$dir/plane-first.json" "$dir/plane-regular.json" --symmetric --training 1470 `
    --output plane-fresh.bin.gz --report plane-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/project_mod_radius.py" `
    --matrix plane-fresh.bin.gz --matrix-report plane-fresh.json --prefix `
    radius-fresh --output radius-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/modular_fibers.py" `
    --output fibers-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/lift_directions.py" `
    --source fibers-fresh.json --output directions-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/recover_anchors.py" `
    --source directions-fresh.json --output anchors-fresh.json
& scripts/serialized_compute.ps1 python -X utf8 -B "$dir/padic_radius.py" --source `
    anchors-fresh.json --level 64 --workers 3 --output padic-fresh.json.gz
~~~

本轮使用本机 Intel Core i5-11400H、16 GiB DDR4，以及已有的 Python 3.14.7、SymPy 1.14.0、mpmath 1.3.0、python-flint 0.9.0、仓库内的 GCC、Boost 和 FLINT。计算作业限制为 6,000,000,000 字节与 50% CPU。全局数学核验主程序运行约 224 秒。此次工作沿用了公开证明代码、精确根数据和仓库既有代数约化，新增部分从四方向参数化开始。
