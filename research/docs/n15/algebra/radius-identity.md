# 2,869 次半径多项式的特征零核验

令 $P$ 表示[整系数文件](../../../n15/algebra/four-directions/reconstruction-p256.json.gz)中按升幂排列的 2,870 个整数。它的次数为 2,869，系数最大位长为 2,931，首项系数为正，全部系数的最大公因数为 1。此文件记录的是从素数幂重建候选的那一步。[四方向消元与全局根来源](four-directions.md)

## 不可约性与实根

对每个保持次数且平方自由的良素数，任何有理因子的次数都必须是模不可约因子次数的子集和。这是 Gauss 引理和模素数后首项保持的直接推论。依次在 65,519、65,521、65,537、65,539、65,543 上分解并精确回乘，五次求交后的可能次数个数依次为 1,482、96、60、12、2。最后剩余 0 和 2,869，故 $P$ 在有理数域不可约。

令 $a=318142930859262836359491181649935005186663866157651845538266/10^{60}$，$b=a+10^{-60}$。有向外包给出 $P(a)>0>P(b)$。将 $P'$ 在中点平移，前 62 个 Taylor 系数各自取球外包。若 $h=(b-a)/2$，未使用部分的绝对值不超过 $2869\,\|P\|_1 2^{2869}h^{62}/(1-h)$，因为中点的模小于 1。包含所有非恒定项的变化量仍小于中点导数绝对值。整个区间上 $P'<0$，故其中恰有一个实根，记为 $\beta$。[整数候选的不可约性与实根证书](../../../n15/algebra/four-directions/integer-candidate-check.json)

## 用五个实坐标表示接触链

使用[四方向推导](four-directions.md)中的单位复方向 $A,B,C,D$，源接点 $E$、左端 $L$ 和右端 $V$。置 $s=(A-1)/(A+1)$，$z=s^2$。对于每个单位量 $Y$ 写

$$ h_Y=\frac{Y-1}{s(Y+1)},\qquad
Y=x_Y+s y_Y,\qquad (x_Y,y_Y)=\left(\frac{1+zh_Y^2}{1-zh_Y^2},\frac{2h_Y}{1-zh_Y^2}\right). $$

目标处 $s$ 纯虚，$z<0$。$h_A=1$。以下依次记 $h_E=e,h_D=d,h_L=l,h_B=b_1$，其中 $b_1$ 与隔离区间右端点无关。半径写作 $t$。虚部等式 $y_E=t y_A$ 给出 $U=t(1-ze^2)-e(1-z)=0$。设

$$ T=ze y_A-x_A-(2x_D-1),\quad x_C=x_B+T/2,\quad \kappa=1-t^2. $$

链的复分式线性矩阵的第一行是 $m=\kappa+t^2 TC$、$n=-t(T+\kappa/B)$，第二行是 $\bar n,\bar m$。设 $m_x=\kappa+t^2 T x_C$，$n_x=-t(T+\kappa x_B)$，$n_y=t\kappa y_B$。将它作用于 $E$ 并取半角坐标，得到

$$ K y_C+J=0,\quad K=t^2T(zle-1),\quad J=l(m_x+n_x-z n_y e)-n_y-(m_x-n_x)e. $$

由 $x_C^2-z y_C^2=1$，可取 $Q=(x_C^2-1)K^2-zJ^2$。在 $K\ne0$ 的目标邻域，它与链有相同的选定分支。

右端的半角 $v=t y_D/(1+t x_D)$。内接点 $W=(w_x,w_y)$ 取 $w_x=1+t(2x_D-1+x_A), w_y=t y_A$，第二坐标仍表示 $s$ 的系数。定义

$$ M=w_x^2-zw_y^2+1,\quad N(h)=M-2w_x+4zw_yh-z(M+2w_x)h^2, $$
$$ J_0=(1+zlv)w_x-z(l+v)w_y+zlv-1,\qquad H=N(l)N(v)-4t^2J_0^2. $$

连接 $L,V$ 的直线方程左端在 $W$ 的值为 $J_0$，其法向量的平方范数为 $(1-zl^2)(1-zv^2)$。又有 $N(h)=(1-zh^2)|W-Y(h)|^2$。由三点外接半径公式，$H=0$ 恰为原几何闭环。程序对 $H,Q$ 各乘已知分母和非零常数，直接构造 $\mathbb Z[t,e,z,d,l,b_1]$ 中的多项式；以下沿用 $H,Q$ 表示清分母后的式子。[整数公式](../../../n15/algebra/four-directions/scalar_incidence.py)

对 $i=e,z,d$ 写 $V_i=H_iQ_l-Q_iH_l$，下标表示固定其余变量的偏导。使用的六个整数式为

$$ U,\quad H,\quad Q,\quad Q_{b_1},\quad V_d,\quad V_z U_e-V_e U_z. $$

原四个方向实参数是半角 $s_A,s_B,s_C,s_D$。$z=-s_A^2,d=s_D/s_A,b_1=s_B/s_A$。在 $s_A\ne0$ 且 $\partial t/\partial s_C\ne0$ 处，$(t,z,d,b_1)$ 是四方向的局部坐标。由 $U_e\ne0,Q_l\ne0$，$U=Q=0$ 局部解出 $e,l$。将其代入 $H$，所得原几何闭环在最优构型处沿固定 $t$ 的三个方向导数均为零。这来自原根已验证的闭环梯度与半径梯度成比例。分别对 $b_1,d,z$ 求导，利用 $H_l\ne0$，恰得最后三个整数式。相关除数和主元可直接在源盒像上验证。[局部唯一性验算](../../../n15/algebra/four-directions/scalar_target.py)

## 用整数特化证明完整的多项式恒等式

实半角坐标用 $x_j=N_j(t)/(c_jP'(t))$ 表示，其中 $c_j$ 是正整数，$N_j$ 是次数小于 2,869 的整系数多项式。重建使用等式

$$ \bar P' x_j\bmod\bar P=\sum_i \frac{\bar P}{M_i}\big(M_i' x_{ij}\bmod M_i\big), $$

这里 $\bar P$ 是模素数幂的首一乘积，$M_i$ 是各支互素的首一半径因子。右边模每个 $M_i$ 等于左边，且次数小于乘积次数。再乘 $P$ 的整数首项并恢复各 $c_j$。以下核验的逻辑只依赖所交付的整数和多项式。

对上述任一整数式 $G$，设 $m$ 为五个坐标的合计次数界，$m_j$ 为各自次数界。把坐标代入，乘 $(P')^m\prod c_j^{m_j}$，得到整系数 $S(t)$。若 $n=\deg P$，则 $\deg S\le\deg_t G+(n-1)m$，且

$$
\|S\|_1\le \|G\|_1\|P'\|_1^m\prod_j
c_j^{m_j}\max\left(1,\frac{\|N_j\|_1}{c_j\|P'\|_1}\right)^{m_j}.
$$

此式按每个单项式使用三角不等式。程序对 $U,H,Q$ 展开统计次数及一范数，对乘积、和与偏导使用整数上界。Sylvester 矩阵的 Hadamard 界给出
$|\operatorname{Res}(P,S)|\le \|P\|_1^{\deg S}\|S\|_1^n$。

选取连续整数 $k_i$，其最大差小于预设整数 $K_0$。从 $|P(k_i)|$ 去除所有不超过 $K_0$ 的素因子，以及整除 $P$ 的首项、$P'(k_i)$ 或任一 $c_j$ 的素因子，留下 $m_i$。在 $\mathbb Z/m_i$ 中直接计算五个坐标并回代六式。如果六式全为零，则每个 $S(k_i)$ 都被 $m_i$ 整除。

**乘积整除引理。** 若 $R=\operatorname{Res}(P,S)\ne0$，则 $\prod_i m_i\mid R$。为证明，固定任意 $q\mid\prod_i m_i$，写 $a_i=v_q(m_i)>0$，仅保留这些下标。因 $q>K_0$，各 $k_i$ 模 $q$ 互异。因 $q$ 不整除首项，$\mathbb Z_q[t]/(P)$ 是秩 $n$ 的自由模。代入各 $k_i$ 给出到 $\bigoplus_i\mathbb Z/q^{a_i}\mathbb Z$ 的满射。它良定义因为 $P(k_i)\equiv0\pmod{q^{a_i}}$，满射由 Lagrange 插值得到，插值分母 $k_i-k_j$ 全是 $q$ 进单位。乘以 $S$ 的像落在核里，因此其余核的阶至少有 $q$ 因子 $q^{\sum_i a_i}$。乘法矩阵的行列式与 $R$ 只差一个 $q$ 进单位，故 $v_q(R)\ge\sum_i a_i$。对各素数使用该结论即得引理，允许不同 $m_i$ 含有共同的素因子。

每个通过的整数点贡献 $\lfloor\log_2 m_i\rfloor$。当总和超过六个结式上界的最大二进制对数，就迫使六个结式均为零。$P$ 不可约，所以 $P$ 整除每个 $S$。根 $\beta$ 处 $P'\ne0$，五个有理坐标满足全部六式。最后将 $P$ 的选定实根进一步以有理端点夹住，把五个有理函数的严格球外包与原构型的唯一根盒比较。[上界与整数点验算](../../../n15/algebra/four-directions/prove_integer_identity.py) · [GMP 验证器](../../../n15/algebra/four-directions/integer_identity.cpp) · [实根连接程序](../../../n15/algebra/four-directions/bind_radius.py)

## 实际核验

五个有理函数分别从 $998244353^{256}$ 和 $998244353^{512}$ 的提升结果恢复，全部分子、整数分母逐项相同。六个清分母方程在有理函数域中的通用公式已经独立展开核对。[五个坐标的全部系数](../../../n15/algebra/four-directions/scalars-p256.json.gz) · [通用公式的符号核验](../../../n15/algebra/four-directions/scalar-identities.json)

从 $k=65,536$ 至 $89,535$ 共 24,000 个整数点逐一通过六式。采用 $K_0=100,000$。各模数贡献的位数之和为 **1,122,055,293**，超过最大的结式上界 **1,120,516,517**。另以独立的前向微分在前三点复核，并故意改动坐标，三个反控制均被拒绝。位数和严格超过上界，故得到六条特征零恒等式。[逐点位数及每式上界](../../../n15/algebra/four-directions/integer-identity.json.gz)

源根像的六维盒半宽为 $10^{-28}$，压缩常数小于 $2.41\times10^{-22}$。对 $P$ 的隔离根使用 13,800 比特有向球算术，将半径和五个坐标分别包在该盒内。前述驻定性说明源根像也是这个六维方程组的零点，唯一性遂证明 $\beta=r_D(15)$。[源盒与主元](../../../n15/algebra/four-directions/scalar-target.json) · [六个坐标的实际夹逼及同根结论](../../../n15/algebra/four-directions/radius-binding.json)
