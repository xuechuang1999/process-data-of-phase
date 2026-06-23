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
```

各模块职责如下：

`vasp_io.py` 只负责从明确的 VASP 文件路径读取原始数据，包括
`OSZICAR`、`CONTCAR` 和 `1-dipol/*/OUTCAR`。

`branch.py` 只负责 Berry phase 分支修正，包括整数晶格矢量平移、逐 image
连接和每个分量连续性检查。

`directions.py` 只负责把 `{001}`、`{110}`、`{111}` 等方向转成单位矢量。

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

其中 \(i\) 是 NEB image 编号，\(U_0=0\)。如果使用 `--energy-source neb`
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

二者单位都是 \(e\cdot\mathrm{\AA}\)。

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

其中 \(\Omega_i\) 的单位是 \(\mathrm{\AA^3}\)。

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

其中 \(n_{i1},n_{i2},n_{i3}\) 是整数。不同 image 需要沿晶轴方向平移不同
整数倍，才能得到物理上连续的极化变化。这个步骤对应
`Pmn2_1-xyz-001-011.xlsx` 中的 `Translate integer unit cell` 行。

当前脚本默认使用 `endpoint-target` 分支选择方式。做法是给定初态和终态的
物理极化目标，并在每个 image 上按反应坐标线性插值得到参考目标：

$$
\mathbf{P}^{\mathrm{target}}_i =
\mathbf{P}^{\mathrm{initial}}+
\frac{i}{N-1}
\left(
\mathbf{P}^{\mathrm{final}}-\mathbf{P}^{\mathrm{initial}}
\right)
$$

然后从 image 00 开始逐 image 展开 Berry 分支。实际选择时，相邻 image 的
polarization 分量连续性优先，\(\mathbf{P}^{\mathrm{target}}_i\) 只作为同等连续
候选之间的弱约束。这样可以避免逐 image 独立贴近目标值时在 image 5、6、7、8
错误跨越一个或多个 polarization quantum。

例如当前 HfO2 路径中，若只用 `center-zero` 并强制右半段单调连接，末态极化
会被展开到约 \(190\ \mu\mathrm{C/cm^2}\)，这是不物理的；若逐 image 独立贴近
`endpoint-target`，image 4 到 image 5 仍会出现约
\(50\ \mu\mathrm{C/cm^2}\) 的跳变。现在默认的连续展开会把相邻 image 的最大
分量跳变控制在约 \(10\ \mu\mathrm{C/cm^2}\) 以内。

这里要区分两个方向：

$$
\hat{\mathbf{b}}
$$

是 Berry 分支展开/单调连接所用的极化翻转方向；而

$$
\hat{\mathbf{e}}
$$

是外电场方向，只在外场能量修正时用于投影。对于当前
`Pca2_1/xyz-000-011` 路径，极化翻转主要沿 z 方向，所以默认设置为：

```bash
--branch-direction 0,0,-1 --branch-coordinates lattice-diagonal
```

即使计算 `{110}` 或 `{111}` 外电场，也不应该让 `{110}` 或 `{111}` 电场方向
反过来决定 Berry 分支；应该先沿极化翻转路径选好分支，再把得到的
\(\Delta\boldsymbol{\mu}_i\) 投影到外电场方向。

默认初态目标为：

$$
\mathbf{P}^{\mathrm{initial}}=(0,0,50)\ \mu\mathrm{C/cm^2}
$$

默认终态目标为：

$$
\mathbf{P}^{\mathrm{final}}=-\mathbf{P}^{\mathrm{initial}}
$$

可通过 `--initial-polarization` 和 `--final-polarization` 修改。例如：

```bash
python compute_coercive_field.py --branch-mode endpoint-target --initial-polarization 0,0,50 --final-polarization 0,0,-50
```

`center-zero` 仍保留为诊断模式。如果需要检查“从 saddle image 附近的 0 极化
分支向两侧延拓”的结果，可以运行：

```bash
python compute_coercive_field.py --branch-mode center-zero --output-dir results_center_zero
```

得到分支修正后的偶极矩后，极化为：

$$
\mathbf{P}_i =
\frac{1600}{\Omega_i}\boldsymbol{\mu}_i
$$

这里 \(\mathbf{P}_i\) 的单位是 \(\mu\mathrm{C/cm^2}\)，
\(\boldsymbol{\mu}_i\) 的单位是 \(e\cdot\mathrm{\AA}\)，
\(\Omega_i\) 的单位是 \(\mathrm{\AA^3}\)。常数 1600 来自单位换算。

相对于初态的极化变化为：

$$
\Delta\mathbf{P}_i =
\frac{1600}{\Omega_i}
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

其中 \(\alpha=x,y,z\)。程序会输出：

$$
\max_{\alpha}
\left|\Delta P_{\alpha,i}^{\mathrm{prev}}\right|
$$

并与 `--continuity-threshold` 给出的阈值比较。默认阈值为
\(80\ \mu\mathrm{C/cm^2}\)。如果某个 image 的任意分量跳变超过阈值，
`branch_correction.csv` 中 `component_continuity_ok` 会是 `False`，
`coercive_field_summary.csv` 的 notes 也会记录对应 image。

## 外电场下的 NEB 能量修正

外电场下采用的自由能修正为：

$$
F_i(\mathbf{E}) =
U_i-\Omega_i\Delta\mathbf{P}_i\cdot\mathbf{E}
$$

这里 \(F_i\) 是相对于 image 00 的外场修正自由能。当前脚本中
\(\Delta\mathbf{P}_i\) 使用 `target_P_*`，也就是从
`--initial-polarization` 到 `--final-polarization` 线性插值得到的切换极化坐标：

$$
\Delta\mathbf{P}_i =
\mathbf{P}^{\mathrm{target}}_i -
\mathbf{P}^{\mathrm{target}}_0
$$

`branch_P_*` 是从 VASP LCALCPOL 结果经 Berry 分支修正得到的诊断量，不再直接
用于外场能量修正。这样 `polarization_profiles.png`、`energy_profiles.png` 和
`field_profiles.png` 使用的是同一条切换极化路径。

实际计算中，电场强度用 \(\mathrm{MV/cm}\)，极化用
\(\mu\mathrm{C/cm^2}\)，体积用 \(\mathrm{\AA^3}\)，能量用 eV，因此脚本使用
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

其中 \(E\) 是电场大小，单位为 \(\mathrm{MV/cm}\)，
\(\hat{\mathbf{e}}\) 是外电场方向的单位矢量。上式等价于更简洁的偶极矩形式：

$$
F_i(E) =
U_i - 0.01 E
\left(\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}\right)
$$

其中 \(E\) 的单位是 \(\mathrm{MV/cm}\)，
\(\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}\) 的单位是
\(e\cdot\mathrm{\AA}\)。

对于每一个给定的电场，外场修正后的 NEB 势垒默认只在中间 image 上取最大值，
也就是排除 image 00 和最后一个 endpoint：

$$
B(E)=\max_{i=1,\ldots,N-1}F_i(E)-F_0(E)
$$

这对应脚本默认参数 `--barrier-mode intermediate`。这样可以避免终态 endpoint
由于数值残余相对能量或几乎为零的投影差控制矫顽场。若要把最后一个 endpoint
也纳入最大值，可以使用：

```bash
python compute_coercive_field.py --barrier-mode full-path
```

若最后一个 endpoint 也纳入，则对应 Excel 中的：

$$
\mathrm{barrier}=\max(F_1,F_2,\ldots,F_N)-F_0
$$

注意这里的最大值不包含 image 00。因此当电场超过临界值后，
\(B(E)\) 可以变成负值；脚本不会把它强行截断为 0。

矫顽场 \(E_c\) 定义为：

$$
B(E_c)=0
$$

因为每个 image 的 \(F_i(E)\) 对电场是线性的，也可以写成：

$$
E_c =
\max_i
\frac{U_i}{
0.01\left(\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}\right)
}
$$

这里只使用
\(\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}>0\) 的 image。给出最大值的
image 就是限制矫顽场的 limiting image。

## 电场方向的设置

Excel 的做法是先给出一个方向整数，例如 \((0,0,-1)\)、\((-1,-1,0)\) 或
\((-1,-1,-1)\)，再用晶格长度归一化。脚本中这个模式叫
`lattice-diagonal`。

设输入方向为：

$$
\mathbf{d}=(d_1,d_2,d_3)
$$

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
--field-coordinates lattice-diagonal
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

## {001}, {110}, {111} 方向下的修正

对于不同外电场方向，零场 NEB 能量 \(U_i\)、分支修正后的偶极矩
\(\boldsymbol{\mu}_i\)、极化变化 \(\Delta\mathbf{P}_i\) 都不需要重新定义。
唯一改变的是把极化变化投影到不同的电场方向：

$$
\Delta\mu^{\parallel}_i =
\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}
$$

随后仍然使用同一个外场修正公式：

$$
F_i(E)=U_i-0.01E\Delta\mu^{\parallel}_i
$$

### {001} 方向

若电场沿负的 \(\{001\}\) 方向：

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
python compute_coercive_field.py --field-direction 0,0,-1 --field-coordinates lattice-diagonal
```

### {110} 方向

若电场沿负的 \(\{110\}\) 方向：

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
python compute_coercive_field.py --field-direction -1,-1,0 --field-coordinates lattice-diagonal --output-dir results_110
```

在 PowerShell 中，由于参数值以 `-` 开头，建议写成：

```bash
python compute_coercive_field.py --field-direction=-1,-1,0 --field-coordinates lattice-diagonal --output-dir results_110
```

### {111} 方向

若电场沿负的 \(\{111\}\) 方向：

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
python compute_coercive_field.py --field-direction -1,-1,-1 --field-coordinates lattice-diagonal --output-dir results_111
```

在 PowerShell 中建议写成：

```bash
python compute_coercive_field.py --field-direction=-1,-1,-1 --field-coordinates lattice-diagonal --output-dir results_111
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

或创建 conda 环境：

```bash
conda env create -f environment.yml
conda activate coercive-field
python compute_coercive_field.py
```

常用变体：

```bash
python compute_coercive_field.py --field-direction 0,0,-1 --output-dir results_001
python compute_coercive_field.py --field-direction=-1,-1,0 --output-dir results_110
python compute_coercive_field.py --field-direction=-1,-1,-1 --output-dir results_111
python compute_coercive_field.py --barrier-mode to-saddle --output-dir results_to_saddle
python compute_coercive_field.py --branch-mode endpoint-target --initial-polarization 0,0,50 --final-polarization 0,0,-50 --output-dir results_endpoint_target
```

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
```

其中：

`coercive_field_summary.csv` 记录每个样品的零场势垒、矫顽场、限制 image、
电场方向、分支展开方向、连续性阈值、初末态极化和 warning notes 等。

`raw_vasp_data.csv` 记录直接从 VASP 文件读取的原始数据。每一行对应一个
sample 的一个 image，包含：

```text
energy_path
structure_path
dipole_outcar_path
energy_eV
rel_energy_eV
lattice_A_*
volume_A3
ion_dipole_eA_*
electronic_dipole_eA_*
raw_dipole_eA_*
raw_P_uC_cm2_*
polarization_quantum_uC_cm2_*
```

这个文件用于检查“原始数据到底从哪个文件读取，以及读到了什么”。

`branch_correction.csv` 记录 Berry 分支修正过程，包含：

```text
branch_n1, branch_n2, branch_n3
branch_P_uC_cm2_*
branch_dipole_eA_*
delta_branch_P_from_prev_uC_cm2_*
component_jump_ok_x/y/z
component_continuity_ok
monotonic_along_branch
```

这个文件用于检查“每个 image 沿每个晶轴方向平移了多少整数倍，以及修正后
每个 polarization 分量是否连续变化”。

`field_projection.csv` 记录外场投影和能量修正相关量，包含：

```text
delta_mu_eA_*
delta_P_uC_cm2_*
delta_mu_parallel_eA
delta_P_parallel_uC_cm2
field_energy_slope_eV_per_MV_cm
required_field_MV_cm
can_be_lowered_by_field
branch_delta_mu_eA_*
branch_delta_P_uC_cm2_*
```

其中 `delta_*` 是实际用于外场能量修正的 `target_P_*` 线性路径投影；
`branch_delta_*` 是 VASP/LCALCPOL 分支修正极化对应的诊断量。这个文件用于检查
“电场方向和 polarization/偶极矩的投影是否正确”。

`polarization_branches.csv` 目前保留为兼容旧输出，内容与
`branch_correction.csv` 相同。

`polarization_profiles.csv` 同时记录两类极化。`target_P_*` 是从
`--initial-polarization` 到 `--final-polarization` 线性插值得到的切换极化坐标；
`branch_P_*` 是从 VASP LCALCPOL 数据经 Berry 分支修正得到的诊断量。
`polarization_profiles.png` 绘制的是 `target_P_switching_uC_cm2`，因此应为一条
随 image 线性变化的直线。

`barrier_vs_field.csv` 记录不同外场下的 \(B(E)\)。

`corrected_profiles.csv` 记录零场和矫顽场下每个 image 的
\(F_i(E)\)。

`field_profiles.png` 是不同电场下的 NEB 曲线族，应该呈现类似 Excel/PDF 中
随外场增大而向终态倾斜的形状。
