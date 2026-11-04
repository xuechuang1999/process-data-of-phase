# 矫顽场计算流程说明

本目录中的脚本 `compute_coercive_field.py` 用来按照论文和
`Pmn2_1-xyz-001-011.xlsx` 的思路，从 VASP 的 NEB 与 LCALCPOL 结果中
计算外电场下的 NEB 自由能曲线，并由势垒随电场降低到 0 的位置得到矫顽场。

## 程序结构

程序已经拆成项目式结构，`compute_coercive_field.py` 只是命令行入口：

```text
compute_coercive_field.py
coercive_field/
  cli.py
  constants.py
  models.py
  vasp_io.py
  branch.py
  directions.py
  calculation.py
  outputs.py
  final_polarization_direction.py
```

各模块职责如下：

`vasp_io.py` 只负责从明确的 VASP 文件路径读取原始数据，包括
`OSZICAR`、`CONTCAR` 和 `1-dipol/*/OUTCAR`。

`branch.py` 只负责 Berry phase 分支修正，包括整数晶格矢量平移、逐 image
连接和每个分量连续性检查。

`directions.py` 只负责把 `{001}`、`{101}`、`{110}`、`{111}` 等晶轴坐标系下的方向转成笛卡尔坐标下的单位矢量。

`calculation.py` 只负责外电场投影、外场自由能修正、势垒和矫顽场计算。

`outputs.py` 统一写出 CSV 和图片，使原始数据和每一步处理结果都可以追踪。

## 需要从 VASP 导出的数据

每一条路径由若干个 NEB image 组成，例如 `00` 到 `08`。对当前默认路径，
未掺杂与掺杂样品分别位于：

```text
Pca2_1/xyz-000-011
Pca2_1/xyz-000-011/a-xyz-000-011
Pca2_1/xyz-000-011/b-xyz-000-011
Pca2_1/xyz-000-011/c-xyz-000-011
```

下面以 `Pca2_1/xyz-000-011` 为例说明文件来源；其它 a/b/c-HZO 子目录完全
相同。对每个 image 需要三类数据。

第一类是零场 NEB 能量。默认从 NEB 路径本身的 image 目录读取：

```text
Pca2_1/xyz-000-011/00/OSZICAR
Pca2_1/xyz-000-011/01/OSZICAR
...
Pca2_1/xyz-000-011/08/OSZICAR
```

脚本读取每个 `OSZICAR` 最后一步的 `E0`，并以 image 00 为参考得到相对能量：

$$
U_i = E_i - E_0
$$

其中 $i$ 是 NEB image 编号，$U_0=0$。如果使用 `--energy-source neb`
或 `--energy-source auto`，也可以从 `neb.dat` 读取已经整理好的相对能量。

注意，能量不从 `1-dipol` 文件夹读取。`1-dipol` 只用于极化相关的静态计算。

第二类是极化相关的偶极矩。需要对每个 NEB image 单独做一次静态计算，
脚本默认读取：

```text
Pca2_1/xyz-000-011/1-dipol/00/OUTCAR
Pca2_1/xyz-000-011/1-dipol/01/OUTCAR
...
Pca2_1/xyz-000-011/1-dipol/08/OUTCAR
```

这些目录中的 `INCAR` 应该是静态极化计算设置。例如当前
`Pca2_1/xyz-000-011/1-dipol/04/INCAR` 中包含：

```text
NSW = 0
IBRION = -1
LCALCPOL=.TRUE.
```

因此 `1-dipol` 中读取的是 LCALCPOL 输出的偶极/极化信息，而不是 NEB 能量。
脚本从 `OUTCAR` 中提取 VASP 打印的两项：

$$
\mathbf{p}^{\mathrm{ion}}_i
$$

$$
\mathbf{p}^{\mathrm{elc}}_i
$$

二者单位都是 $e\cdot\mathrm{\AA}$。

第三类是晶格信息。脚本从 NEB 路径本身每个 image 的 `CONTCAR` 读取晶格
矢量和体积：

```text
Pca2_1/xyz-000-011/00/CONTCAR
Pca2_1/xyz-000-011/01/CONTCAR
...
Pca2_1/xyz-000-011/08/CONTCAR
```

不是从 `1-dipol` 文件夹读取晶格，也不是优先用 `POSCAR`。`CONTCAR` 是 NEB
弛豫后的每个 image 结构；只有在缺少 `CONTCAR` 时，脚本才会退回到同一
image 目录下的 `POSCAR`。

$$
\mathbf{a}_i,\mathbf{b}_i,\mathbf{c}_i
$$

$$
\Omega_i = |\det(\mathbf{a}_i,\mathbf{b}_i,\mathbf{c}_i)|
$$

其中 $\Omega_i$ 的单位是 $\mathrm{\AA^3}$。

## 从 VASP 数据得到极化分支

首先把离子偶极矩和电子偶极矩相加，得到每个 image 的原始总偶极矩：

$$
\boldsymbol{\mu}^{\mathrm{raw}}_i =
\mathbf{p}^{\mathrm{ion}}_i+\mathbf{p}^{\mathrm{elc}}_i
$$

Berry phase 极化有量子不唯一性，所以需要给总偶极矩加上整数个晶格矢量，
选择一条连续的极化分支。Excel 中的 `Translate integer unit cell` 行就是
这些整数。脚本中对应为：

$$
\boldsymbol{\mu}_i =
\boldsymbol{\mu}^{\mathrm{raw}}_i+
n_{i1}\mathbf{a}_i+n_{i2}\mathbf{b}_i+n_{i3}\mathbf{c}_i
$$

其中 $n_{i1},n_{i2},n_{i3}$ 是整数。不同 image 需要沿晶轴方向平移不同
整数倍，才能得到物理上连续的极化变化。这个步骤对应
`Pmn2_1-xyz-001-011.xlsx` 中的 `Translate integer unit cell` 行。

当前脚本默认使用 `center-zero` 分支选择方式。它不需要用户指定逐 image 的
target polarization 来决定分支，而是分三步完成：

第一步，先把相邻 image 放在同一个连续 Berry 分支上。对每个 image 枚举整数
分支：

$$
\mathbf{P}^{\mathrm{candidate}}_i =
\mathbf{P}^{\mathrm{raw}}_i+
n_{i1}\mathbf{Q}_{a,i}+
n_{i2}\mathbf{Q}_{b,i}+
n_{i3}\mathbf{Q}_{c,i}
$$

其中 $\mathbf{Q}_{a,i},\mathbf{Q}_{b,i},\mathbf{Q}_{c,i}$ 是该 image 沿
$\mathbf{a},\mathbf{b},\mathbf{c}$ 三个方向的 polarization quantum。程序选择
使相邻 image 分量跳变最小的整数 $n_{i1},n_{i2},n_{i3}$，从而让整条路径位于
同一个连续 branch 上。

第二步，在已经连续的 branch 上，对所有 image 施加同一个整数平移
$(k_1,k_2,k_3)$：

$$
\mathbf{P}^{\mathrm{branch}}_i =
\mathbf{P}^{\mathrm{continuous}}_i+
k_1\mathbf{Q}_{a,i}+
k_2\mathbf{Q}_{b,i}+
k_3\mathbf{Q}_{c,i}
$$

这个全局整数平移由中间 image 决定，使中间 image 的 $|\mathbf{P}^{\mathrm{branch}}|$ 最小，也就是最靠近 0。它不要求中间 image 严格等于 0，因为 polarization quantum 是离散的。

第三步，若 `--final-polarization auto` 根据路径旋转矩阵识别到末态不是 180 度反向，
程序会用初态极化和旋转矩阵得到末态极化参考值，再只把中间 image 之后到末态的
右半段整体平移到对应 Berry branch。这一步不是把所有 image 贴到线性 target
path 上，而是只解除右端 branch 简并。

这里要区分两个概念。

第一，Berry 分支矫正本身是三维矢量矫正。默认 `center-zero` 模式不会指定
逐 image 的 target polarization，也不会只沿某一个方向平移；它会同时考虑
$\mathbf{a},\mathbf{b},\mathbf{c}$ 三个方向的 polarization quantum，让
$\mathbf{P}^{\mathrm{branch}}_i$ 在三个分量上都连续。

第二，`--field-direction` 给出的是外电场方向：

$$
\hat{\mathbf{e}}
$$

它只在外场能量修正时用于点乘投影。对于当前
`Pca2_1/xyz-000-011` 路径，铁电极化方向是 $\{001\}$。若研究沿初始铁电极化
反向的外场，默认设置为：

```bash
--field-direction 0,0,-1 --field-coordinates lattice
```

`branch_correction.csv` 里的 `branch_projection_uC_cm2` 和
`monotonic_along_branch` 诊断列直接使用外电场方向投影。当前不再提供单独的
branch 方向参数。

## Pca2_1 与 Pmn2_1 的方向

这两类相的铁电极化方向和 NEB 末态方向不能混为一谈。

- `Pca2_1` 相的铁电极化方向是 $\{001\}$。但 NEB 计算的末态可以不是 180 度
  反向翻转；例如可以从原来的 $\{001\}$ 方向转到 $\{010\}$ 方向。
- `Pmn2_1` 相的铁电极化方向是 $\{101\}$。NEB 末态也可以是非 180 度转换；
  例如可以从原来的 $\{101\}$ 方向转到 $[-110]$ 方向。

因此，默认分支矫正不能写成“指定一个初末态 target polarization，然后强制贴近
这条 target 线”。正确做法是先让 LCALCPOL 给出的矢量极化在 Berry branch 上
连续，再整体平移 polarization quantum，使中间 image 最接近 0。

如果只是为了输出 `target_P_*` 诊断参考列，可以按实际路径给
`--initial-polarization` 和 `--final-polarization`。例如 `Pca2_1` 的
非 180 度路径可以让初态与 $\{001\}$ 成比例，末态与 $\{010\}$ 成比例；
`Pmn2_1` 的非 180 度路径可以让初态与 $\{101\}$ 成比例，末态与 $[-110]$
成比例。这里的三个数应当是实际极化矢量分量，单位为
$\mu\mathrm{C/cm^2}$，Miller 指数只说明方向，不能直接替代极化大小。

外场方向也要单独指定。例如沿 `Pca2_1` 初始铁电极化反向的场可写成：

```bash
python compute_coercive_field.py --field-direction 0,0,-1 --field-coordinates lattice
```

沿 `Pmn2_1` 的 $\{101\}$ 反向场可写成：

```bash
python compute_coercive_field.py --field-direction=-1,0,-1 --field-coordinates lattice
```

若要研究非 180 度末态方向上的外场，例如 `Pca2_1` 的 $\{010\}$ 或
`Pmn2_1` 的 $[-110]$，则分别把外场方向改为 `0,1,0` 或 `-1,1,0`。
这只是改变 $\hat{\mathbf{e}}$；分支矫正后的
$\Delta\boldsymbol{\mu}_i$ 仍然来自同一条矢量极化路径。

`--initial-polarization` 和 `--final-polarization` 仍保留，用来生成
`target_P_*` 诊断列。例如默认诊断参考为：

$$
\mathbf{P}^{\mathrm{initial}}=(0,0,50)\ \mu\mathrm{C/cm^2}
$$

并且：

```bash
python compute_coercive_field.py --initial-polarization 0,0,50 --final-polarization 0,0,-50
```

注意，这条 `target_P_*` 参考线不参与默认外场能量修正；外场能量修正使用
Berry 分支矫正后的 $\mathbf{P}^{\mathrm{branch}}_i$。

得到分支修正后的偶极矩后，极化为：

$$
\mathbf{P}_i =
\frac{C_{\mu P}}{\Omega_i}\boldsymbol{\mu}_i
$$

这里 $\mathbf{P}_i$ 的单位是 $\mu\mathrm{C/cm^2}$，
$\boldsymbol{\mu}_i$ 的单位是 $e\cdot\mathrm{\AA}$，
$\Omega_i$ 的单位是 $\mathrm{\AA^3}$。当前代码中
`UC_CM2_PER_EANG_A3` 取精确换算值：

$$
C_{\mu P}=1602.176634
$$

偶极矩除以体积的单位为：

$$
\frac{e\cdot\mathrm{\AA}}{\mathrm{\AA^3}}
=
\frac{e}{\mathrm{\AA^2}}
$$

换成 SI 单位：

$$
1\frac{e\cdot\mathrm{\AA}}{\mathrm{\AA^3}}
=
\frac{
1.602176634\times10^{-19}\ \mathrm{C}
\times 10^{-10}\ \mathrm{m}
}{
10^{-30}\ \mathrm{m^3}
}
=
16.02176634\ \mathrm{C/m^2}
$$

又因为：

$$
1\ \mathrm{C/m^2}
=
100\ \mu\mathrm{C/cm^2}
$$

所以：

$$
1\frac{e\cdot\mathrm{\AA}}{\mathrm{\AA^3}}
=
1602.176634\ \mu\mathrm{C/cm^2}
$$

因此脚本使用：

$$
\mathbf{P}\ [\mu\mathrm{C/cm^2}]
=
\frac{1602.176634}{\Omega\ [\mathrm{\AA^3}]}
\boldsymbol{\mu}\ [e\cdot\mathrm{\AA}]
$$

反过来，从目标极化换算目标偶极矩变化时使用：

$$
\Delta\boldsymbol{\mu}\ [e\cdot\mathrm{\AA}]
=
\frac{
\Delta\mathbf{P}\ [\mu\mathrm{C/cm^2}]
\Omega\ [\mathrm{\AA^3}]
}{1602.176634}
$$

如果为了和旧 Excel 手算近似对照，也常把 $C_{\mu P}$ 写成 1600；当前源码使用
精确值 1602.176634。

相对于初态的极化变化为：

$$
\Delta\mathbf{P}_i =
\frac{C_{\mu P}}{\Omega_i}
\left(\boldsymbol{\mu}_i-\boldsymbol{\mu}_0\right)
$$

等价地，也可以直接使用偶极矩差：

$$
\Delta\boldsymbol{\mu}_i =
\boldsymbol{\mu}_i-\boldsymbol{\mu}_0
$$

这些中间量会写入：

```text
results/branch_correction.csv
```

其中包含分支修正后的 polarization、translated dipole、整数分支编号、相邻
image 的每个分量变化量，以及连续性检查结果。

相邻 image 的分量连续性检查使用：

$$
\Delta P_{\alpha,i}^{\mathrm{prev}} =
P_{\alpha,i}-P_{\alpha,i-1}
$$

其中 $\alpha=x,y,z$。程序会输出：

$$
\max_{\alpha}
\left|\Delta P_{\alpha,i}^{\mathrm{prev}}\right|
$$

并与 `--continuity-threshold` 给出的阈值比较。默认 `auto` 阈值由
当前路径的 polarization quantum 自动得到：

$$
P_{\mathrm{th}} =
\frac{
\max_i\max_{\mathbf{m}\in\{-1,0,1\}^3\setminus\{\mathbf{0}\}}
\left|\mathbf{m}\cdot\mathbf{Q}_i\right|
}{2N_{\mathrm{image}}}
$$

其中 $\mathbf{Q}_i=(\mathbf{Q}_{a,i},\mathbf{Q}_{b,i},\mathbf{Q}_{c,i})$。
这样对 `101` 这类斜向翻转会考虑 $\mathbf{Q}_a+\mathbf{Q}_c$ 这样的组合极化量子。
如果某个 image 的任意分量跳变超过阈值，
`branch_correction.csv` 中 `component_continuity_ok` 会是 `False`，
`coercive_field_summary.csv` 的 notes 也会记录对应 image。

## 外电场下的 NEB 能量修正

外电场下采用的自由能修正为：

$$
F_i(\mathbf{E}) =
U_i-\Omega_i\Delta\mathbf{P}_i\cdot\mathbf{E}
$$

这里 $F_i$ 是相对于 image 00 的外场修正自由能。当前脚本中
$\Delta\mathbf{P}_i$ 使用 Berry 分支修正后的矢量极化 `branch_P_*`：

$$
\Delta\mathbf{P}_i =
\mathbf{P}^{\mathrm{branch}}_i -
\mathbf{P}^{\mathrm{branch}}_0
$$

也可以等价地使用分支修正后的偶极矩差：

$$
\Delta\boldsymbol{\mu}_i =
\boldsymbol{\mu}^{\mathrm{branch}}_i -
\boldsymbol{\mu}^{\mathrm{branch}}_0
$$

因此外场能量修正使用的是一个真正的矢量点乘：

$$
\Delta\mathbf{P}_i\cdot\mathbf{E}
$$

`target_P_*` 仍会输出，但只是从 `--initial-polarization` 到
`--final-polarization` 的线性参考线，用于诊断和画图，不参与默认外场能量修正。

这个公式里的符号、物理量和程序数据来源如下。

- $i$：NEB image 编号，对应 `00`、`01` 等 image 目录，也对应各个 CSV 中的
  `image` 列。
- $E_i$：第 $i$ 个 image 的零场能量，默认来自该 image 的 `OSZICAR` 最后一步
  `E0`；若使用 `--energy-source neb` 或可用的 `auto`，相对能量来自 `neb.dat`。
- $U_i=E_i-E_0$：相对于 image 00 的零场 NEB 能量，对应
  `raw_vasp_data.csv` 和 `field_projection.csv` 中的 `rel_energy_eV`。
- $\Omega_i$：第 $i$ 个 image 的晶胞体积，来自该 image 的 `CONTCAR`
  或缺省回退的 `POSCAR`，对应 `raw_vasp_data.csv` 中的 `volume_A3`。
- $\mathbf{P}^{\mathrm{branch}}_i$：VASP LCALCPOL 极化经过 Berry 分支修正后的
  矢量极化，对应 `branch_correction.csv` 中的 `branch_P_uC_cm2_*`。
- $\Delta\mathbf{P}_i=\mathbf{P}^{\mathrm{branch}}_i-\mathbf{P}^{\mathrm{branch}}_0$：
  相对于 image 00 的分支修正极化变化，对应 `field_projection.csv` 中的
  `delta_P_uC_cm2_*`。
- $\boldsymbol{\mu}^{\mathrm{branch}}_i$：与 `branch_P_*` 对应的分支修正偶极矩，
  对应 `branch_dipole_eA_*`。
- $\Delta\boldsymbol{\mu}_i$：相对于 image 00 的分支修正偶极矩变化，对应
  `field_projection.csv` 中的 `delta_mu_eA_*`。
- $\mathbf{E}$：外电场矢量。程序命令行没有直接输入完整矢量
  $\mathbf{E}$，而是输入电场大小 $E$ 和方向 $\hat{\mathbf{e}}$，即
  $\mathbf{E}=E\hat{\mathbf{e}}$。$E$ 的单位是 $\mathrm{MV/cm}$，
  $\hat{\mathbf{e}}$ 由 `--field-direction` 和 `--field-coordinates` 得到。
- $F_i(E)$：外场修正后的相对自由能，对应 `corrected_profiles.csv` 中的
  `corrected_rel_energy_eV`。

注意，矢量公式 $-\Omega_i\Delta\mathbf{P}_i\cdot\mathbf{E}$ 和后面的投影写法
不是两次修正。因为 $\mathbf{E}=E\hat{\mathbf{e}}$，所以：

$$
\Delta\mathbf{P}_i\cdot\mathbf{E}
=
E\left(\Delta\mathbf{P}_i\cdot\hat{\mathbf{e}}\right)
$$

括号中的 $\Delta\mathbf{P}_i\cdot\hat{\mathbf{e}}$ 只是把极化变化取到外电场
方向上的分量；垂直于外电场的极化分量对 $-\mathbf{P}\cdot\mathbf{E}$ 没有贡献。
程序把这个分量写入 `delta_P_parallel_uC_cm2`，这样同一条分支修正极化路径可以方便地
换成 `{001}`、`{110}`、`{111}` 等不同外场方向。

实际计算中，电场强度用 $\mathrm{MV/cm}$，极化用
$\mu\mathrm{C/cm^2}$，体积用 $\mathrm{\AA^3}$，能量用 eV，因此脚本使用
下面的单位换算形式：

$$
F_i(E) =
U_i-
\frac{
\left(\Delta\mathbf{P}_i\cdot\hat{\mathbf{e}}\right)
E
\times 10^8
\times \Omega_i
\times 10^{-30}
}{
100\times 1.602176634\times 10^{-19}
}
$$

其中 $E$ 是电场大小，单位为 $\mathrm{MV/cm}$，
$\hat{\mathbf{e}}$ 是外电场方向的单位矢量。上式等价于更简洁的偶极矩形式：

$$
F_i(E) =
U_i - 0.01 E
\left(\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}\right)
$$

其中 $E$ 的单位是 $\mathrm{MV/cm}$，
$\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}$ 的单位是
$e\cdot\mathrm{\AA}$。

默认矫顽场判据为逐段下坡判据，对应参数：

```bash
--coercive-criterion downhill --barrier-mode full-path
```

也就是要求路径上每一步后一个 image 的外场自由能都不高于前一个 image：

$$
F_i(E)-F_{i-1}(E)\le 0
$$

程序输出和 `barrier_vs_field.csv` 中的 `barrier_eV` 在这个模式下表示当前
电场下路径上最大的相邻上坡量：

$$
D(E)=\max_{i=1,\ldots,N}\left[F_i(E)-F_{i-1}(E)\right]
$$

矫顽场 $E_c$ 定义为：

$$
D(E_c)=0
$$

因为每个 image 的 $F_i(E)$ 对电场是线性的，所以逐段下坡判据对应：

$$
E_c =
\max_i
\frac{U_i-U_{i-1}}{
0.01\left[
\left(\Delta\boldsymbol{\mu}_i-\Delta\boldsymbol{\mu}_{i-1}\right)
\cdot\hat{\mathbf{e}}
\right]
}
$$

这里只使用
$\left(\Delta\boldsymbol{\mu}_i-\Delta\boldsymbol{\mu}_{i-1}\right)\cdot\hat{\mathbf{e}}>0$ 且 $U_i-U_{i-1}>0$ 的相邻段。给出最大值的
相邻段就是限制矫顽场的 limiting constraint，例如 `1->2`。

这里 $N$ 是最后一个 image 的编号；$D(E)$ 是给定外场下最大的相邻上坡能量，
对应默认 `downhill` 模式下 `barrier_vs_field.csv` 的 `barrier_eV`；
$E_c$ 是让最大相邻上坡刚好降到 0 的电场强度，对应
`coercive_field_summary.csv` 中的 `coercive_field_MV_cm`。

旧的整体势垒判据仍可使用：

```bash
python compute_coercive_field.py --coercive-criterion barrier --barrier-mode intermediate
```

它对应：

$$
B(E)=\max_iF_i(E)-F_0(E)
$$

这里 $B(E)$ 是整体势垒，也就是选定 image 范围内最高的外场自由能相对
image 00 的高度。使用 `--coercive-criterion barrier` 时，
`barrier_vs_field.csv` 的 `barrier_eV` 表示 $B(E)$。

以及：

$$
E_c =
\max_i
\frac{U_i-U_0}{
0.01\left[
\left(\Delta\boldsymbol{\mu}_i-\Delta\boldsymbol{\mu}_0\right)
\cdot\hat{\mathbf{e}}
\right]
}
$$

## 电场方向的设置

Excel 的做法是先给出一个方向整数，例如 $(0,0,-1)$、$(-1,-1,0)$ 或
$(-1,-1,-1)$，再用晶格长度归一化。脚本中这个模式叫
`lattice`。

设输入方向为：

$$
\mathbf{d}=(d_1,d_2,d_3)
$$

这里 $\mathbf{d}$ 是命令行 `--field-direction` 给出的方向标签或方向矢量，
例如 `{001}` 方向常写成 $(0,0,-1)$。它本身还不是单位矢量。

对于近似正交晶胞，Excel 使用的归一化方式为：

$$
L_{\mathbf{d}}=
\sqrt{
(a_xd_1)^2+(b_yd_2)^2+(c_zd_3)^2
}
$$

$$
\hat{\mathbf{e}}=
\frac{
(a_xd_1,\ b_yd_2,\ c_zd_3)
}{
L_{\mathbf{d}}
}
$$

这对应脚本参数：

```bash
--field-coordinates lattice
```

如果想用完整晶格矢量做分数坐标到笛卡尔坐标的转换，则使用：

$$
\hat{\mathbf{e}}=
\frac{
d_1\mathbf{a}+d_2\mathbf{b}+d_3\mathbf{c}
}{
|d_1\mathbf{a}+d_2\mathbf{b}+d_3\mathbf{c}|
}
$$

对应参数：

```bash
--field-coordinates lattice
```

如果输入的方向已经是笛卡尔方向，则使用：

```bash
--field-coordinates cartesian
```

## 不同外电场方向下的修正

对于不同外电场方向，零场 NEB 能量 $U_i$、分支修正后的矢量极化
$\mathbf{P}^{\mathrm{branch}}_i$，以及由它换算得到的
$\Delta\boldsymbol{\mu}_i$ 都不需要重新定义。
唯一改变的是把偶极矩变化投影到不同的外电场方向：

$$
\Delta\mu^{\parallel}_i =
\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}
$$

随后仍然使用同一个外场修正公式：

$$
F_i(E)=U_i-0.01E\Delta\mu^{\parallel}_i
$$

这里 $\Delta\mu^{\parallel}_i$ 是偶极矩变化沿外电场方向的标量分量，
对应 `field_projection.csv` 中的 `delta_mu_parallel_eA`。因此“改变方向”
不是重新定义 NEB 路径或极化路径，而是改变 $\hat{\mathbf{e}}$，从而改变
$\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}$。

### 方向

若电场沿负的 $\{001\}$ 方向：

$$
\mathbf{d}_{001}=(0,0,-1)
$$

则：

$$
\hat{\mathbf{e}}_{001}=
\frac{(0,0,-c_z)}{|c_z|}
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction 0,0,-1 --field-coordinates lattice
```

这通常用于 `Pca2_1` 相沿 $\{001\}$ 铁电极化轴的反向外场。

### 方向

`Pmn2_1` 相的铁电极化方向是 $\{101\}$。若电场沿负的 $\{101\}$ 方向：

$$
\mathbf{d}_{101}=(-1,0,-1)
$$

则：

$$
\hat{\mathbf{e}}_{101}=
\frac{(-a_x,0,-c_z)}{\sqrt{a_x^2+c_z^2}}
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction=-1,0,-1 --field-coordinates lattice --output-dir results_101
```

### 方向

`Pca2_1` 的某些 NEB 路径末态可能转到 $\{010\}$ 方向。若要研究沿这个末态
方向的外场，可以使用：

$$
\mathbf{d}_{010}=(0,1,0)
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction 0,1,0 --field-coordinates lattice --output-dir results_010
```

### [-110] 方向

`Pmn2_1` 的某些 NEB 路径末态可能转到 $[-110]$ 方向。按常见晶向写法，
$[-110]$ 表示 $(-1,1,0)$，也就是只让第一个指数取负。若要研究沿这个末态
方向的外场，可以使用：

$$
\mathbf{d}_{[-110]}=(-1,1,0)
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction=-1,1,0 --field-coordinates lattice --output-dir results_m110
```

### 负的 [110] 方向

若电场沿负的 $[110]$ 方向，即 $-[110]=[-1,-1,0]$：

$$
\mathbf{d}_{110}=(-1,-1,0)
$$

则：

$$
\hat{\mathbf{e}}_{110}=
\frac{(-a_x,-b_y,0)}{\sqrt{a_x^2+b_y^2}}
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction -1,-1,0 --field-coordinates lattice --output-dir results_110
```

在 PowerShell 中，由于参数值以 `-` 开头，建议写成：

```bash
python compute_coercive_field.py --field-direction=-1,-1,0 --field-coordinates lattice --output-dir results_110
```

这里和 $[-110]$ 不是同一个方向；前者是 $(-1,-1,0)$，后者是 $(-1,1,0)$。
如果某个 Excel 文件的 `Direction` 行写的是 `-1,-1,0`，就按这一行的三个整数
直接传给 `--field-direction`。

### 方向

若电场沿负的 $\{111\}$ 方向：

$$
\mathbf{d}_{111}=(-1,-1,-1)
$$

则：

$$
\hat{\mathbf{e}}_{111}=
\frac{(-a_x,-b_y,-c_z)}{\sqrt{a_x^2+b_y^2+c_z^2}}
$$

运行命令为：

```bash
python compute_coercive_field.py --field-direction -1,-1,-1 --field-coordinates lattice --output-dir results_111
```

在 PowerShell 中建议写成：

```bash
python compute_coercive_field.py --field-direction=-1,-1,-1 --field-coordinates lattice --output-dir results_111
```

如果要完全复现某个 Excel 文件，请优先查看该表中 `Direction` 行的三个整数，
并把它们直接传给 `--field-direction`。

## 输入目录格式

默认样品路径为：

```text
Pca2_1/xyz-000-011
Pca2_1/xyz-000-011/a-xyz-000-011
Pca2_1/xyz-000-011/b-xyz-000-011
Pca2_1/xyz-000-011/c-xyz-000-011
```

每条路径应包含：

```text
00/CONTCAR ... 08/CONTCAR
00/OSZICAR ... 08/OSZICAR
1-dipol/00/OUTCAR ... 1-dipol/08/OUTCAR
1-dipol/00/INCAR ... 1-dipol/08/INCAR
```

其中 `00/CONTCAR ... 08/CONTCAR` 和 `00/OSZICAR ... 08/OSZICAR`
来自 NEB 路径本身；`1-dipol/00/OUTCAR ... 1-dipol/08/OUTCAR`
来自对每个 NEB image 结构进行的 LCALCPOL 静态计算。

`OSZICAR` 是默认能量来源。`neb.dat` 可以通过下面的方式启用：

```bash
python compute_coercive_field.py --energy-source neb
```

或：

```bash
python compute_coercive_field.py --energy-source auto
```

## 运行方式

直接使用当前 Python 环境：

```bash
python compute_coercive_field.py
```

或创建 `conda` 环境：

```bash
conda env create -f environment.yml
conda activate coercive-field
python compute_coercive_field.py
```

常用变体：

```bash
python compute_coercive_field.py --field-direction 0,0,-1 --output-dir results_001
python compute_coercive_field.py --field-direction 0,1,0 --output-dir results_010
python compute_coercive_field.py --field-direction=-1,0,-1 --output-dir results_101
python compute_coercive_field.py --field-direction=-1,1,0 --output-dir results_m110
python compute_coercive_field.py --field-direction=-1,-1,0 --output-dir results_110
python compute_coercive_field.py --field-direction=-1,-1,-1 --output-dir results_111
python compute_coercive_field.py --barrier-mode to-saddle --output-dir results_to_saddle
python compute_coercive_field.py --coercive-criterion barrier --barrier-mode intermediate --output-dir results_barrier
python compute_coercive_field.py --initial-polarization 0,0,50 --final-polarization 0,0,-50 --output-dir results_manual_polarization
```

批量脚本 `run_all_fe2_paths.sh` 会默认扫描 `001`、`010`、`101`、`m110`
和 `111` 五个外场方向，其中 `010` 对应 `Pca2_1` 常见非 180 度末态，
`m110` 对应 `Pmn2_1` 的 $[-110]$ 末态方向 `(-1,1,0)`。

## 用户可控制参数

下面列出 `compute_coercive_field.py` 的所有命令行参数。方向或矢量参数均写成
`x,y,z`，不要加空格。若在 PowerShell 中输入负数开头的方向，建议使用等号形式，
例如 `--field-direction=-1,0,-1`。

### 输入与样品

| 参数                | 默认值                 | 可选值或格式                   | 含义                                                                                                                                                 |
| ------------------- | ---------------------- | ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--root`          | `Pca2_1/xyz-000-011` | 路径                           | 主 NEB 路径目录。若不指定`--sample`，程序会自动读取该目录以及同级的 `a-*`、`b-*`、`c-*` 子目录。                                             |
| `--sample`        | 无                     | `LABEL=PATH`，可重复         | 手动指定要计算的样品。只要使用了任意`--sample`，默认的 undoped/a/b/c 自动样品列表就不再使用。                                                      |
| `--energy-source` | `oszicar`            | `oszicar`、`neb`、`auto` | 零场能量来源。`oszicar` 从每个 image 的 `OSZICAR` 读取；`neb` 从 `neb.dat` 读取；`auto` 优先尝试 `neb.dat`，不可用时回退到 `OSZICAR`。 |

示例：

```bash
python compute_coercive_field.py --root="Pca2_1\xyz-000-011"
python compute_coercive_field.py --sample HfO2=Pca2_1\xyz-000-011 --sample a-HZO=Pca2_1\xyz-000-011\a-xyz-000-011
```

### 极化与 Berry Branch

| 参数                       | 默认值   | 可选值或格式                                          | 含义                                                                                                                                                                                                                                |
| -------------------------- | -------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--initial-polarization` | `auto` | `auto` 或 `Px,Py,Pz`，单位 $\mu\mathrm{C/cm^2}$ | `auto` 会先从路径识别相结构初态方向：`Pca2_1 -> 001`，`Pmn2_1 -> 101`，方向按晶轴坐标解释；再从 VASP LCALCPOL 的 image 00 等价 branch 中选择与该方向一致的初态极化。手动传入时仍作为线性`target_P_*` 诊断参考线的初态极化。 |
| `--final-polarization`   | `auto` | `auto`、`opposite-initial` 或 `Px,Py,Pz`        | 线性`target_P_*` 诊断参考线的末态极化。`auto` 会根据路径旋转关系从初态极化推断；`opposite-initial` 表示取 `--initial-polarization` 的相反数。                                                                               |
| `--continuity-threshold` | `auto` | `auto` 或浮点数，单位$\mu\mathrm{C/cm^2}$         | 判断相邻 image 三个极化分量是否连续的阈值；`auto` 使用 $\max|\mathbf{m}\cdot\mathbf{Q}|/(2N_{\mathrm{image}})$，超过则在输出 notes 和`component_continuity_ok` 中标记。                                                       |

常见设置：

```bash
python compute_coercive_field.py --initial-polarization auto --final-polarization auto
python compute_coercive_field.py --final-polarization opposite-initial
python compute_coercive_field.py --final-polarization 0,60,0
```

### 外电场与矫顽场

| 参数                     | 默认值        | 可选值或格式                                   | 含义                                                                                                                                |
| ------------------------ | ------------- | ---------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `--field-direction`    | `0,0,-1`    | `reverse-initial` 或 `dx,dy,dz`            | 外电场方向$\hat{\mathbf{e}}$。外场修正使用分支修正后的 $\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}$。                        |
| `--field-coordinates`  | `lattice`   | `lattice`、`cartesian`                     | 解释`--field-direction` 的坐标方式。默认按完整晶格基矢解释晶轴坐标，例如 `-1,0,0` 表示沿 $-\mathbf{a}$。                      |
| `--barrier-mode`       | `full-path` | `to-saddle`、`intermediate`、`full-path` | 参与势垒或下坡判据的 image 范围。`to-saddle` 只到零场 saddle；`intermediate` 排除两端；`full-path` 排除 image 00 但包含末态。 |
| `--coercive-criterion` | `downhill`  | `downhill`、`barrier`                      | 矫顽场判据。`downhill` 要求每个相邻段都不再上坡；`barrier` 使用整体势垒降为 0 的判据。                                          |

示例：

```bash
python compute_coercive_field.py --field-direction 0,0,-1 --field-coordinates lattice
python compute_coercive_field.py --field-direction=-1,0,-1 --field-coordinates lattice
python compute_coercive_field.py --coercive-criterion barrier --barrier-mode intermediate
```

### 输出与电场扫描

| 参数             | 默认值      | 可选值或格式                   | 含义                                                                  |
| ---------------- | ----------- | ------------------------------ | --------------------------------------------------------------------- |
| `--field-step` | `0.1`     | 浮点数，单位$\mathrm{MV/cm}$ | `barrier_vs_field.csv` 和 `barrier_vs_field.png` 的电场扫描步长。 |
| `--field-max`  | 自动        | 浮点数，单位$\mathrm{MV/cm}$ | 电场扫描上限。不指定时，程序根据计算得到的最大有限$E_c$ 自动选择。  |
| `--output-dir` | `results` | 路径                           | 输出 CSV 和 PNG 的目录。                                              |

示例：

```bash
python compute_coercive_field.py --field-step 0.05 --field-max 40 --output-dir results_fine_scan
```

### 批量脚本参数

`run_all_fe2_paths.sh` 会对一个 phase 目录下所有形如 `xyz-101-000` 的路径批量运行
`compute_coercive_field.py`。它自身可控制的参数如下。

| 参数                    | 默认值                                      | 含义                                                                                              |
| ----------------------- | ------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| `--phase-root DIR`    | 脚本所在目录的上一级                        | 包含所有 FE2 路径文件夹的 phase 目录。                                                            |
| `--code-dir DIR`      | 脚本所在目录                                | 包含`compute_coercive_field.py` 的代码目录。                                                    |
| `--output-root DIR`   | `<phase-root>/coercive_field_all_results` | 批量输出根目录。                                                                                  |
| `--python CMD`        | `python` 或环境变量 `$PYTHON_BIN`       | 调用 Python 的命令。                                                                              |
| `--extra-args "ARGS"` | 空                                          | 额外传给`compute_coercive_field.py` 的参数，例如 `--extra-args "--field-direction=-1,0,-1"`。 |
| `--dry-run`           | 关闭                                        | 只打印将要运行的命令，不实际计算。                                                                |
| `-h`、`--help`      | 无                                          | 显示帮助。                                                                                        |

## 当前默认 example 的实际处理结果

本节给出当前源码按默认命令：

```bash
python compute_coercive_field.py
```

对默认 example 重新分析得到的结果。默认参数等价于：

```text
root = Pca2_1/xyz-000-011
samples = undoped-HfO2, a-HZO, b-HZO, c-HZO
energy_source = oszicar
field_direction = 0,0,-1
field_coordinates = lattice
barrier_mode = full-path
coercive_criterion = downhill
initial_polarization = auto
final_polarization = auto
continuity_threshold = auto
UC_CM2_PER_EANG_A3 = 1602.176634
```

当前默认路径会读取四个 sample：

| label            | sample path                          | image |
| ---------------- | ------------------------------------ | ----: |
| `undoped-HfO2` | `Pca2_1/xyz-000-011`               | 00-08 |
| `a-HZO`        | `Pca2_1/xyz-000-011/a-xyz-000-011` | 00-08 |
| `b-HZO`        | `Pca2_1/xyz-000-011/b-xyz-000-011` | 00-08 |
| `c-HZO`        | `Pca2_1/xyz-000-011/c-xyz-000-011` | 00-08 |

每个 image 读取并处理的变量链如下：

| 变量                                       | 来源                                          | 处理方式                                                                                                                         | 输出位置                                                               |
| ------------------------------------------ | --------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| $E_i$                                    | `00/OSZICAR ... 08/OSZICAR` 最后一个 `E0` | 以 image 00 为参考，$U_i=E_i-E_0$                                                                                              | `raw_vasp_data.csv`: `rel_energy_eV`                               |
| $\mathbf{a}_i,\mathbf{b}_i,\mathbf{c}_i$ | `00/CONTCAR ... 08/CONTCAR`                 | 组成晶格矩阵                                                                                                                     | `raw_vasp_data.csv`: `lattice_A_*`                                 |
| $\Omega_i$                               | 同一个`CONTCAR`                             | $\Omega_i=\lvert\det(\mathbf{a}_i,\mathbf{b}_i,\mathbf{c}_i)\rvert$                                                            | `raw_vasp_data.csv`: `volume_A3`                                   |
| $\mathbf{p}^{ion}_i,\mathbf{p}^{elc}_i$  | `1-dipol/00/OUTCAR ... 1-dipol/08/OUTCAR`   | VASP LCALCPOL 偶极矩，单位$e\cdot\mathrm{\AA}$                                                                                 | `raw_vasp_data.csv`: `ion_dipole_eA_*`, `electronic_dipole_eA_*` |
| $\boldsymbol{\mu}^{raw}_i$               | 上一行两项                                    | $\boldsymbol{\mu}^{raw}_i=\mathbf{p}^{ion}_i+\mathbf{p}^{elc}_i$                                                               | `raw_vasp_data.csv`: `raw_dipole_eA_*`                             |
| $\mathbf{P}^{raw}_i$                     | $\boldsymbol{\mu}^{raw}_i,\Omega_i$         | $\mathbf{P}^{raw}_i=C_{\mu P}\boldsymbol{\mu}^{raw}_i/\Omega_i$                                                                | `raw_vasp_data.csv`: `raw_P_uC_cm2_*`                              |
| $\mathbf{Q}_{a,b,c,i}$                   | 晶格和体积                                    | $\mathbf{Q}_{a,i}=C_{\mu P}\mathbf{a}_i/\Omega_i$ 等                                                                           | `raw_vasp_data.csv`: `polarization_quantum_uC_cm2_*`               |
| $n_{i1},n_{i2},n_{i3}$                   | Berry 分支搜索                                | 选整数，使$\mathbf{P}^{branch}_i=\mathbf{P}^{raw}_i+n_{i1}\mathbf{Q}_{a,i}+n_{i2}\mathbf{Q}_{b,i}+n_{i3}\mathbf{Q}_{c,i}$ 连续 | `branch_correction.csv`: `branch_n1/2/3`                           |
| $\mathbf{P}^{branch}_i$                  | 原始极化和整数分支                            | 分支修正后的 LCALCPOL 极化，用于外场修正和分支诊断                                                                               | `branch_correction.csv`: `branch_P_uC_cm2_*`                       |
| $\mathbf{P}^{target}_i$                  | 初末态目标极化诊断参考                        | 从$(0,0,50)$ 到 $(0,0,-50)$ 线性插值，不参与默认外场修正                                                                     | `polarization_profiles.csv`: `target_P_uC_cm2_*`                   |
| $\Delta\mathbf{P}_i$                     | 分支修正极化变化                              | $\Delta\mathbf{P}_i=\mathbf{P}^{branch}_i-\mathbf{P}^{branch}_0$                                                               | `field_projection.csv`: `delta_P_uC_cm2_*`                         |
| $\Delta\boldsymbol{\mu}_i$               | 分支修正偶极矩变化                            | $\Delta\boldsymbol{\mu}_i=\boldsymbol{\mu}^{branch}_i-\boldsymbol{\mu}^{branch}_0$                                             | `field_projection.csv`: `delta_mu_eA_*`                            |
| $\Delta\mu_i^\parallel$                  | 外场方向                                      | $\Delta\mu_i^\parallel=\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}$，当前 $\hat{\mathbf{e}}=(0,0,-1)$                      | `field_projection.csv`: `delta_mu_parallel_eA`                     |
| $F_i(E)$                                 | $U_i,\Delta\mu_i^\parallel$                 | $F_i(E)=U_i-0.01E\Delta\mu_i^\parallel$                                                                                        | `corrected_profiles.csv`: `corrected_rel_energy_eV`                |
| $E_c$                                    | 相邻下坡约束                                  | 对每段$i-1\to i$，$E_i^{req}=(U_i-U_{i-1})/[0.01(\Delta\mu_i^\parallel-\Delta\mu_{i-1}^\parallel)]$，取最大值                | `coercive_field_summary.csv`: `coercive_field_MV_cm`               |

四个 sample 的最终结果为：

| label            | saddle image | zero-field barrier (eV) | $E_c$ (MV/cm) | limiting constraint | endpoint$U$ (eV) | endpoint$\Delta\mu^\parallel$ ($e\cdot\mathrm{\AA}$) | notes                                                                             |
| ---------------- | -----------: | ----------------------: | --------------: | ------------------- | -----------------: | -------------------------------------------------------: | --------------------------------------------------------------------------------- |
| `undoped-HfO2` |            4 |                0.330640 |       19.553880 | 1->2                |           0.000000 |                                                 5.084166 | auto endpoint branch direction 0,1,0                                              |
| `a-HZO`        |            4 |                0.251390 |       13.628965 | 1->2                |           0.001500 |                                                 5.098636 | auto endpoint branch direction 0,1,0                                              |
| `b-HZO`        |            4 |                0.248530 |       13.160084 | 1->2                |           0.000500 |                                                 5.144821 | auto endpoint branch direction 0,1,0                                              |
| `c-HZO`        |            4 |                0.251440 |       13.391692 | 1->2                |           0.004210 |                                                 5.120286 | `neb.dat` 与 `OSZICAR` 最大差 0.0697 eV；auto endpoint branch direction 0,1,0 |

四个 sample 的 Berry 分支诊断端点为：

| label            | $\mathbf{P}^{branch}_{00}$ ($\mu\mathrm{C/cm^2}$) | $\mathbf{P}^{branch}_{08}$ ($\mu\mathrm{C/cm^2}$) | 最大相邻分量跳变 ($\mu\mathrm{C/cm^2}$) | 连续性 |
| ---------------- | ----------------------------------------------------- | ----------------------------------------------------- | ----------------------------------------: | ------ |
| `undoped-HfO2` | (-0.001, 0.000, 50.057)                               | (-0.001, 59.769, -10.064)                             |                                    59.931 | True   |
| `a-HZO`        | (0.003, -0.006, 50.485)                               | (-0.014, 59.304, -9.045)                              |                                    59.419 | True   |
| `b-HZO`        | (0.011, 3.812, 50.411)                                | (0.001, 55.626, -9.655)                               |                                    52.010 | True   |
| `c-HZO`        | (3.180, -0.017, 50.557)                               | (-3.051, 59.258, -9.233)                              |                                    59.434 | True   |

注意：上表中的 `branch_P_*` 是 LCALCPOL/Berry 分支矫正后的矢量极化。当前默认
外场能量修正直接使用 `branch_P_*` 对应的 `branch_dipole_*`，再与外电场方向点乘。
`target_P_*` 只作为线性参考线输出。

下面把 `undoped-HfO2` 这个 sample 逐 image 展开，展示每一步从读取数据到外场约束的处理结果。
向量均按 $(x,y,z)$ 写出；`req_E` 是当前 `downhill` 判据下该相邻段单独需要的电场。

| image | $U_i$ (eV) | $\Omega_i$ ($\mathrm{\AA^3}$) | $\boldsymbol{\mu}^{raw}_i$ ($e\cdot\mathrm{\AA}$) | $\mathbf{P}^{raw}_i$ ($\mu\mathrm{C/cm^2}$) | branch$n$  | $\mathbf{P}^{branch}_i$ ($\mu\mathrm{C/cm^2}$) | $\Delta\mathbf{P}^{branch}_i$ ($\mu\mathrm{C/cm^2}$) | $\Delta\mu^\parallel_i$ ($e\cdot\mathrm{\AA}$) | $\Delta E_{step}$ (eV) | req_E (MV/cm) |
| ----: | -----------: | --------------------------------: | ----------------------------------------------------- | ----------------------------------------------- | ------------ | -------------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------: | -----------------------: | ------------: |
|    00 |     0.000000 |                           135.489 | (-0.000, -60.653, 44.906)                             | (-0.001, -717.230, 531.024)                     | (0, 12, -8)  | (-0.001, 0.000, 50.057)                            | (0.000, 0.000, 0.000)                                    |                                              0.000 |                      nan |           nan |
|    01 |     0.059760 |                           136.451 | (-0.000, -60.869, 44.247)                             | (-0.000, -714.703, 519.535)                     | (0, 12, -8)  | (0.000, -0.002, 40.824)                            | (0.001, -0.002, -9.232)                                  |                                              0.756 |                    0.060 |         7.903 |
|    02 |     0.193150 |                           136.550 | (0.000, -61.003, 43.576)                              | (0.000, -715.768, 511.284)                      | (0, 12, -8)  | (-0.000, -0.004, 32.791)                           | (0.001, -0.004, -17.266)                                 |                                              1.438 |                    0.133 |        19.554 |
|    03 |     0.305450 |                           136.012 | (0.000, -61.043, 42.831)                              | (0.001, -719.063, 504.533)                      | (0, 12, -8)  | (0.000, -0.011, 24.639)                            | (0.001, -0.012, -25.418)                                 |                                              2.141 |                    0.112 |        15.973 |
|    04 |     0.330640 |                           135.350 | (0.000, -61.016, 42.107)                              | (0.001, -722.264, 498.438)                      | (0, 12, -8)  | (-0.000, -0.018, 16.800)                           | (0.001, -0.018, -33.257)                                 |                                              2.814 |                    0.025 |         3.746 |
|    05 |     0.305810 |                           136.015 | (0.000, 61.040, 42.824)                               | (0.001, 719.016, 504.438)                       | (0, -11, -9) | (0.000, 59.913, -35.434)                           | (0.001, 59.912, -85.491)                                 |                                              7.241 |                   -0.025 |         0.000 |
|    06 |     0.194090 |                           136.550 | (0.000, 61.002, 43.569)                               | (0.000, 715.756, 511.203)                       | (0, -11, -9) | (0.000, 59.645, -27.104)                           | (0.001, 59.644, -77.161)                                 |                                              6.543 |                   -0.112 |         0.000 |
|    07 |     0.059670 |                           136.452 | (-0.000, 60.868, 44.249)                              | (-0.001, 714.696, 519.556)                      | (0, -11, -9) | (-0.001, 59.559, -18.991)                          | (-0.001, 59.559, -69.048)                                |                                              5.851 |                   -0.134 |         0.000 |
|    08 |     0.000000 |                           135.489 | (-0.000, 60.653, 44.906)                              | (-0.001, 717.231, 531.024)                      | (0, -11, -9) | (-0.001, 59.769, -10.064)                          | (0.000, 59.769, -60.121)                                 |                                              5.084 |                   -0.060 |         0.000 |

这个表可以按下面方式读：

- `raw_dipole` 和 `raw_P` 来自 `1-dipol/*/OUTCAR` 与对应 image 的体积换算，尚未处理 Berry 分支不唯一性。
- `branch_n` 是为了让 `branch_P` 沿 image 连续而添加的整数 polarization quantum。
- `Delta P_branch` 是分支修正后的矢量极化相对于 image 00 的变化。
- `Delta_mu_parallel` 是 `branch_dipole` 相对于 image 00 的变化再投影到外场方向 $(0,0,-1)$ 的结果。
- `req_E` 最大的一段是 `1->2`，因此 `undoped-HfO2` 的限制约束是 `1->2`，最终
  $E_c=19.553880\ \mathrm{MV/cm}$。

## 输出文件

默认输出目录为 `results/`：

```text
coercive_field_summary.csv
raw_vasp_data.csv
branch_correction.csv
field_projection.csv
barrier_vs_field.csv
corrected_profiles.csv
polarization_branches.csv
polarization_profiles.csv
barrier_vs_field.png
energy_profiles.png
field_profiles.png
polarization_profiles.png
branch_polarization_profiles.png
```

列名约定：

- 后缀 `_x/_y/_z` 表示笛卡尔三维矢量的 $x,y,z$ 分量。
- `lattice_A_a_x` 这类列表示晶格矢量矩阵；例如 `lattice_A_b_z` 是
  $\mathbf{b}$ 晶格矢量的 $z$ 分量，单位为 $\mathrm{\AA}$。
- `polarization_quantum_uC_cm2_a_x` 这类列表示沿 $\mathbf{a}$、$\mathbf{b}$、
  $\mathbf{c}$ 平移一个晶胞对应的 polarization quantum 矢量分量。

### `coercive_field_summary.csv`

每一行对应一个 sample，是最终结果总表。

- `label`：样品标签，例如 `undoped-HfO2`、`a-HZO`。
- `path`：该 sample 的 NEB 路径目录。
- `notes`：读取数据、分支连续性、无法降低的约束等 warning 信息。
- `branch_mode`：Berry 分支选择模式。
- `barrier_mode`：参与势垒或下坡判据的 image 范围。
- `coercive_criterion`：矫顽场判据，`downhill` 或 `barrier`。
- `continuity_threshold_uC_cm2`：判断相邻 image 极化分量跳变是否过大的阈值。
- `saddle_image`：零场相对能量最高的 image 编号。
- `zero_field_barrier_eV`：零场最高相对能量，单位 eV/cell。
- `coercive_field_MV_cm`：计算得到的矫顽场 $E_c$，单位 $\mathrm{MV/cm}$。
- `limiting_image`：决定 $E_c$ 的 image；`downhill` 模式下是相邻段的后一个 image。
- `limiting_constraint`：决定 $E_c$ 的约束；`downhill` 模式下写成 `i-1->i`。
- `endpoint_rel_energy_eV`：最后一个 image 相对 image 00 的零场能量。
- `endpoint_projected_dipole_eA`：最后一个 image 的分支修正偶极矩变化在外场方向上的投影。
- `field_dir_*`：外电场单位方向 $\hat{\mathbf{e}}$。
- `branch_dir_*`：用于 `branch_projection_*` 和 `monotonic_along_branch` 的诊断投影方向；
  默认 `center-zero` 分支选择不由这个方向决定。
- `P_initial_uC_cm2_*`、`P_final_uC_cm2_*`：VASP/LCALCPOL 分支修正后的初态和末态极化。
- `target_P_initial_uC_cm2_*`、`target_P_final_uC_cm2_*`：线性目标极化参考线的初态和末态，仅用于诊断。

### `raw_vasp_data.csv`

每一行对应一个 sample 的一个 image，记录从 VASP 文件直接读取或直接换算得到的数据。

- `label`、`image`：样品标签和 NEB image 编号。
- `sample_path`：sample 的 NEB 路径目录。
- `energy_path`：该 image 的 `OSZICAR` 路径。
- `structure_path`：该 image 的结构文件路径，优先 `CONTCAR`，缺失时用 `POSCAR`。
- `dipole_outcar_path`：该 image 的 LCALCPOL 静态计算 `OUTCAR` 路径。
- `energy_eV`：当前写出的所选能量来源数值。代码在读入时已经把能量转成
  image 00 参考下的相对能量，因此通常与 `rel_energy_eV` 相同。
- `rel_energy_eV`：相对 image 00 的零场能量 $U_i$，单位 eV。
- `volume_A3`：晶胞体积 $\Omega_i$，单位 $\mathrm{\AA^3}$。
- `lattice_A_*`：晶格矢量 $\mathbf{a}_i,\mathbf{b}_i,\mathbf{c}_i$，单位 $\mathrm{\AA}$。
- `polarization_quantum_uC_cm2_*`：各晶格方向的 polarization quantum，单位 $\mu\mathrm{C/cm^2}$。
- `ion_dipole_eA_*`：LCALCPOL 输出的离子偶极矩，单位 $e\cdot\mathrm{\AA}$。
- `electronic_dipole_eA_*`：LCALCPOL 输出的电子偶极矩，单位 $e\cdot\mathrm{\AA}$。
- `raw_dipole_eA_*`：原始总偶极矩，`ion_dipole_eA_* + electronic_dipole_eA_*`。
- `raw_P_uC_cm2_*`：由原始总偶极矩和体积换算得到的原始极化。

### `branch_correction.csv`

每一行对应一个 sample 的一个 image，记录 Berry 分支修正过程。这个文件用于检查
每个 image 沿晶轴方向平移了多少个 polarization quantum，以及修正后是否连续。

- `label`、`image`：样品标签和 NEB image 编号。
- `branch_n1`、`branch_n2`、`branch_n3`：沿 $\mathbf{a}$、$\mathbf{b}$、
  $\mathbf{c}$ 方向加入的整数分支编号。
- `branch_projection_uC_cm2`：分支修正后极化在 `branch_dir_*` 诊断方向上的投影。
- `delta_branch_projection_from_prev_uC_cm2`：上述投影相对前一个 image 的变化。
- `max_abs_component_jump_from_prev_uC_cm2`：相邻 image 之间三个极化分量跳变绝对值的最大值。
- `component_jump_ok_x/y/z`：对应分量的相邻跳变是否小于连续性阈值。
- `monotonic_along_branch`：沿 `branch_dir_*` 诊断方向的投影是否满足单调性检查。
- `component_continuity_ok`：三个分量是否都通过连续性检查。
- `branch_P_uC_cm2_*`：Berry 分支修正后的极化。
- `branch_dipole_eA_*`：与 `branch_P_uC_cm2_*` 对应的分支修正后偶极矩。
- `delta_branch_P_from_prev_uC_cm2_*`：分支修正后极化相对前一个 image 的变化。

`polarization_branches.csv` 目前保留为兼容旧输出，内容与 `branch_correction.csv` 相同。

### `field_projection.csv`

每一行对应一个 sample 的一个 image，记录外场方向投影和矫顽场约束。这里的
`delta_*` 是实际用于外场能量修正的 Berry 分支修正路径；`target_delta_*` 是
线性目标路径的诊断量。

- `label`、`image`：样品标签和 NEB image 编号。
- `rel_energy_eV`：零场相对能量 $U_i$。
- `coercive_criterion`：该行按 `downhill` 或 `barrier` 判据解释。
- `delta_mu_eA_*`：分支修正偶极矩变化 $\Delta\boldsymbol{\mu}_i=
  \boldsymbol{\mu}^{branch}_i-\boldsymbol{\mu}^{branch}_0$。
- `delta_P_uC_cm2_*`：分支修正极化变化 $\Delta\mathbf{P}_i=
  \mathbf{P}^{branch}_i-\mathbf{P}^{branch}_0$。
- `delta_mu_parallel_eA`：$\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}$。
- `delta_P_parallel_uC_cm2`：$\Delta\mathbf{P}_i\cdot\hat{\mathbf{e}}$。
- `field_energy_slope_eV_per_MV_cm`：$F_i(E)$ 对电场强度的斜率，
  等于 $-0.01\,\Delta\mu_i^\parallel$。
- `constraint_from_image`、`constraint_to_image`：当前判据检查的约束两端；
  `downhill` 模式为相邻段，`barrier` 模式为 image 00 到当前 image。
- `constraint_delta_energy_eV`：该约束两端的零场能量差。
- `constraint_delta_mu_parallel_eA`：该约束两端的投影偶极矩差。
- `constraint_delta_P_parallel_uC_cm2`：该约束两端的投影极化差。
- `constraint_field_energy_slope_eV_per_MV_cm`：该约束能量差随电场变化的斜率。
- `required_field_MV_cm`：把该约束刚好降到 0 所需的电场。
- `can_be_lowered_by_field`：外场是否能降低该约束；若为 `False`，说明投影方向与降低能垒不匹配。
- `branch_delta_mu_eA_*`：与 `delta_mu_eA_*` 相同，保留这个列名方便识别其来源。
- `branch_delta_P_uC_cm2_*`：与 `delta_P_uC_cm2_*` 相同，保留这个列名方便识别其来源。
- `target_delta_mu_eA_*`：线性目标路径对应的偶极矩变化，仅用于诊断。
- `target_delta_P_uC_cm2_*`：线性目标路径对应的极化变化，仅用于诊断。

### `polarization_profiles.csv`

每一行对应一个 sample 的一个 image，方便画极化路径和比较目标路径与 VASP 分支修正路径。

- `label`、`image`：样品标签和 NEB image 编号。
- `P_parallel_field_uC_cm2`：分支修正极化 `branch_P_*` 在外场方向上的投影。
- `P_parallel_branch_uC_cm2`：分支修正极化在 `branch_dir_*` 诊断方向上的投影。
- `target_P_switching_uC_cm2`：目标极化在切换轴方向上的投影，仅作为 CSV 诊断列保留。
- `target_P_parallel_field_uC_cm2`：目标极化在外场方向上的投影。
- `target_P_parallel_branch_uC_cm2`：目标极化在 `branch_dir_*` 诊断方向上的投影。
- `target_P_uC_cm2_*`：由初末态目标极化线性插值得到的目标极化。
- `branch_P_uC_cm2_*`：VASP/LCALCPOL 数据经 Berry 分支修正后的极化。

### `barrier_vs_field.csv`

每一行对应一个 sample 的一个扫描电场值。

- `label`：样品标签。
- `field_MV_cm`：外电场强度 $E$，单位 $\mathrm{MV/cm}$。
- `barrier_eV`：当前判据值。`downhill` 模式下是最大相邻上坡量 $D(E)$；
  `barrier` 模式下是整体势垒 $B(E)$。
- `path_value_eV`：与 `barrier_eV` 相同，保留这个列名是为了强调它代表当前路径判据值。
- `coercive_criterion`：该行使用的判据。

### `corrected_profiles.csv`

每一行对应一个 sample、一个电场值、一个 image。

- `label`：样品标签。
- `field_MV_cm`：外电场强度。默认只写出 0 和计算得到的 $E_c$。
- `image`：NEB image 编号。
- `corrected_rel_energy_eV`：外场修正后的相对自由能 $F_i(E)$。

### 图片文件

- `barrier_vs_field.png`：不同外场下的 `barrier_eV` 曲线。
- `energy_profiles.png`：零场 NEB 能量曲线。
- `field_profiles.png`：多个外场下的 $F_i(E)$ 曲线族。
- `polarization_profiles.png`：VASP/LCALCPOL 分支修正极化 `branch_P_*` 在外场方向上的投影。
- `branch_polarization_profiles.png`：VASP/LCALCPOL 分支修正极化的 $x/y/z$ 三个分量，用来检查 branch 连续性。
