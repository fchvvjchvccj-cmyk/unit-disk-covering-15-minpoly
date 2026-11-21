# 用 15 个等半径闭圆盘覆盖单位闭圆盘

设 $P(t)=\sum_{j=0}^{2869}a_jt^j$。[系数文件](../../n15/solution/coefficients.json)逐项列出 $a_0,a_1,\ldots,a_{2869}$，每行一个十进制整数。最优半径是 $P$ 在下列有理开区间中的唯一实根 $\alpha$：

```math
\left(
\frac{318142930859262836359491181649935005186663866157651845538266}{10^{60}},
\frac{318142930859262836359491181649935005186663866157651845538267}{10^{60}}
\right).
```

## 最小多项式与构型的连接

$P$ 本原且首项为正。五个良素数上各自分解后，逐个验证模因子的不可约性并回乘，所有可能有理因子次数的交集为 $\{0,2869\}$。$P$ 在左端为正、右端为负，其导数在整段严格为负。以上两项分别证明有理不可约性和单变量根隔离。[全部模因子与根数检查](../../n15/algebra/four-directions/integer-candidate-check.json)

原构型在实半角坐标 $(t,e,z,d,l,b_1)$ 中满足六个整数方程。原构型的像严格落入给定的有理盒，盒上通过了严格自映和压缩检验。另一方面，把 $P$ 的指定实根代入五个明确的有理函数，也得到该盒内的一点。对各函数清分母后，用整数特化和结式的整数上界证明六个余式全部被 $P$ 整除。这一步证明该点精确满足六式。盒中只有一个根，故它的半径坐标 $\alpha$ 就是原构型的半径。[六式的推导与结式证明](algebra/radius-identity.md) · [恒等式计算](../../n15/algebra/four-directions/integer-identity.json.gz) · [实根与源根的连接](../../n15/algebra/four-directions/radius-binding.json)

## 十五个有序圆心

令 $\xi$ 为[整数方程组](../../n15/proofs/global-final/replay/core/system15.py)在[有理盒](../../n15/proofs/global-final/replay/core/root15.json)中的唯一根。盒中心是逐项的 `xnum/Qx`，每一维半宽 $10^{-80}$。从零编号，$\xi_{76}=\alpha^2$。所交付构型依次为

```math
c_{i+1}=(\xi_{20+2i},\xi_{21+2i}),\qquad 0\le i<15.
```

[有序圆心数据](../../n15/solution/ordered-centers.json)逐行给出对应的坐标索引和有理盒中心分子。预条件后的压缩映射严格把盒映入内部，唯一性适用于整组 148 个坐标，因而也唯一指定这组有序圆心。[原根核验程序](../../n15/proofs/global-final/replay/core/verify_root15.py) · [148 元根至消元输入的连接](../../n15/proofs/global-final/axis-chain-20261118.json)

## 连续覆盖与全局最优性

对于上述构型，39 个有向三角形拼成以 11 个边界接点为顶点的拓扑圆盘，绕数检验给出内部的覆盖；逐面用公共点和双圆三角形引理验证。外侧余下的 11 个圆弓通过弧段引理验证，合起来覆盖整个闭单位盘。因此 $R_D(C_*)\le\alpha$。[覆盖证明及引理](proofs/global-final/source/PROOF-reading.zh.md) · [覆盖核验程序](../../n15/proofs/global-final/replay/core/verify_upper15.py)

若存在半径严格小于 $\alpha$ 的覆盖，正松量允许一般位置扰动。取极小子覆盖后，Voronoi 对偶归约把反例送入五族辅助杆图，此归约处理重合圆心、冗余圆盘及退化情形。独立组合计数以轨道质量和根面递推验证了五族的完备性。在 11,950,884 个组合类型中，度量不等式排除 11,907,870 个。剩余的 43,013 个非候选类型和一个候选类型经过穷尽分支核验，未决叶节点为零。于是对任意 $C\in(\mathbb R^2)^{15}$ 均有 $R_D(C)\ge\alpha$，从而 $R_D(C_*)=r_D(15)=\alpha$。[反例归约](proofs/global-final/reduction.zh.md) · [全局下界证明](proofs/global-final/global-proof.zh.md) · [完整重放结果](../../n15/proofs/global-final/master-20261118.json)

## 复核与来源

### 公开资料核对

从公开仓库可直接查看[完整证明](../../../PROOF.zh.md)、[整数系数](../../n15/solution/coefficients.json)、[有序圆心](../../n15/solution/ordered-centers.json)、[实根区间与协议指纹](../../n15/solution/solution.json)及[已保存的复核结果](../../n15/solution/replay-verified.json)。查看这些材料、核对系数和摘要，不等同于从头运行数学证明。

[v1.0.0 发布页](https://github.com/fchvvjchvccj-cmyk/unit-disk-covering-15-minpoly/releases/tag/v1.0.0)上的自定义离线 ZIP 与 GitHub 自动生成的源码 ZIP 不同。完整离线 ZIP 的 SHA-256 为 `656f61628e2868b1e0003d92c5d716168127d3cd4d20c6da7387cf7d4d4996da`；解压后可按包内逐文件 SHA-256 清单核对内容。这些摘要用于核验文件完整性，不能替代数学重放。尚未在仅有该 ZIP 的环境中验证全量重放。

### 全量重放：原研究环境

以下命令记录原研究环境的复算方式，并非全新克隆即可运行的说明。若从公开仓库根目录开始，先进入 `research`，再检查原输入提交是否可用：

~~~powershell
Set-Location research
git cat-file -t ecadc602d1f93c74d528d222dae24178ac255e6c
~~~

只有检查返回 `commit` 且所需 Windows 工具链已经就绪时，再从 `research` 运行以下命令。输出目录必须尚未创建；若示例目录已存在，请换一个新目录。

~~~powershell
& .\scripts\serialized_compute.ps1 python -X utf8 -B n15/solution/verify.py --output-directory .tools/n15-solution-fresh
~~~

目前公开 GitHub 仓库无法解析上述固定提交。它包含五份压缩枚举输入与固定版本的上游程序；重放器从 Git 读取所需文件，在临时目录展开并校验摘要。因此，修正工作目录并不能补齐这些输入。

原 Windows 环境使用 Python 3.14.7、SymPy 1.14.0、mpmath 1.3.0、python-flint 0.9.0，以及 GCC、Boost、FLINT 和 GMP。重放器核对固定编译器及 Boost 头文件摘要；整数验证器使用本机的 `.tools/msolve-paged/chain/ucrt64` 工具链。原计算限制为 6 GB 内存和 50% CPU。完整复核会重跑全局证明、不可约性、实根隔离、整数恒等式及实根连接，并逐字节比对导出的整数列表、圆心与协议指纹。[固定输入和编译器说明](../../n15/proofs/global-final/replay_archived_inputs.py) · [本次消元及环境记录](algebra/four-directions.md)

有序根数据、连续覆盖和组合全局下界沿用固定版本的公开证明。四方向代数约化、半径多项式重建、不可约性、根隔离及特征零连接在本次研究中完成。协议指纹按两个隔离端点分别运行官方脚本，两次舍入和哈希相同。[固定版上游证明包](https://github.com/VonEquinox/DiskCoveringSolve/tree/089584973f41e6374ce110f6e5b5ee456df19352/cover15/proof_bundle) · [交付清单与指纹](../../n15/solution/solution.json)

## 实际运行记录

从本次目标启动到数学核验完成并导出指纹，自然经过约 **4 小时 8 分钟**，包含等待与检验。任务系统在导出后整理报告时显示累计 **4 小时 8 分 21 秒**。本次属于使用公开来源的本地研究，仓库合同的时间限额为 12 小时。固定题目元数据的参考指纹状态为 `not_published`，公开基准成绩未登记。

| 从目标启动起 | 实际结果 |
|---|---|
| 约 9 分钟 | 全局证明重放完成，未决叶节点为零。 |
| 约 32 分钟 | 从 148 元全局根到 25 元代数输入的独立连接完成。 |
| 约 57 分钟 | 直接整数消元遇到进程内存限制后，改用二次代数，取得包含 644,413 项的完整三变量整数表。 |
| 约 1 小时 51 分钟 | 算出只含半径的有限域结式，并对其重数一层逐因子回代。 |
| 约 3 小时 7 分钟 | 从素数幂恢复 2,870 个整数系数，随后证实不可约性及隔离区间的实根数。 |
| 约 3 小时 56 分钟 | 两种精度恢复出完全相同的五个有理函数，并验证选定实根的坐标落入源根盒。 |
| 约 4 小时 8 分钟 | 24,000 个整数特化完成最终特征零证明，随后导出系数、圆心和官方协议指纹。 |

本机为 Intel Core i5-11400H、16 GiB DDR4。核心计算使用已有工具，包括 MSYS2 GCC 16.2.0 和带 FLINT 3.6.0 的 python-flint 0.9.0。全局重放主程序耗时约 224 秒；末级整数恒等式验算耗时约 526 秒。独立核验包括重新编译的公开证明、两条方程路线的交叉检查、独立微分、破坏坐标的反控制，以及对官方指纹的独立编码核对。第三方人工审查和公开参考哈希比对未计入这些检查。模型运行费用与精确推理配置未记录。
