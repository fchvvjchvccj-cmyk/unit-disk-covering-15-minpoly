# 十五个等半径圆盘覆盖单位圆盘的最小多项式

用 15 个等半径闭圆盘覆盖半径为 1 的闭圆盘，允许圆心任意放置。最小公共半径约为 **0.318142930859262836359491181649935**。它在有理数域上的次数为 **2,869**。本仓库给出半径本身的最小多项式、全部 **2,870 个整数系数**、一个严格隔离区间，以及指定实根等于几何最优值的计算机辅助证明。

2026-11-21 修订的中文论文按摘要、引言、预备知识、主要结果、结论和参考文献组织正文。主定理位于第 2 页，后续各节给出连续覆盖上界、任意圆心配置的全局下界，以及最优半径与多项式隔离实根的对应。[阅读八页论文《十五个全等圆盘覆盖单位圆盘的最优半径》](十五个全等圆盘覆盖单位圆盘的最优半径.pdf)

[阅读完整证明](PROOF.zh.md) · [查看七页 PDF](十五圆盘最优覆盖证明.pdf) · [下载完整离线资料包](https://github.com/fchvvjchvccj-cmyk/unit-disk-covering-15-minpoly/releases/download/v1.0.0/unit-disk-covering-15-minpoly-v1.0.0.zip)

## 证明与数据

| 材料 | 内容 |
| --- | --- |
| [整系数数组](research/n15/solution/coefficients.json) | 从常数项到最高次项排列。本原、首项为正；除以首项系数即得首一最小多项式。 |
| [有序圆心](research/n15/solution/ordered-centers.json) | 由多项式系统及有理隔离盒精确指定 15 个圆心。 |
| [实根、隔离端点和协议指纹](research/n15/solution/solution.json) | 指定区间宽度为 10⁻⁶⁰。 |
| [同根证明](research/docs/n15/algebra/radius-identity.md) | 将隔离实根提升到经过唯一性认证的六维几何根盒。 |
| [复核结果](research/n15/solution/replay-verified.json) | 11,950,884 个组合类型全部核对，未决叶节点为零；24,000 个整数特化点验证代数恒等式。 |
| [证明附录](attachments/n15-proof-appendix-20261118.zip) | 文稿所附的精确方程、根盒、五个有理函数、整数证书与推导。 |

## 来源与复核

几何构型、覆盖整个闭单位盘的上界以及排除所有更小半径覆盖的全局下界沿用 VonEquinox 的固定版本，并已在原项目重放。四方向代数约化、2,869 次半径多项式的重建、不可约性、实根隔离及同根认证见本仓库的完整证明。[上游几何证明包](https://github.com/VonEquinox/DiskCoveringSolve/tree/089584973f41e6374ce110f6e5b5ee456df19352/cover15/proof_bundle)

下载并解压 **v1.0.0** 的完整 ZIP，打开 `打开文稿.html` 即可离线阅读。包内含原始附件、上游目录全部 75 个文件、复核代码、历史输入索引和逐文件 SHA-256 清单。复核程序及其所需的固定 Git 输入和 Windows 工具链见[复核入口](research/n15/solution/verify.py)与[运行记录](research/docs/n15/solution.zh.md)。

**完整离线 ZIP 的 SHA-256**

```text
656f61628e2868b1e0003d92c5d716168127d3cd4d20c6da7387cf7d4d4996da
```

**按 diskcover-minpoly-v1 协议计算的答案 SHA-256**

```text
27d58d1f801408165315e66d2d9fcc823b76aa5531581432c84c31dbbe163335
```
