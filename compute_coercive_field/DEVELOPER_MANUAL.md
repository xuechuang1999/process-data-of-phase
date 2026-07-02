# 开发者手册：矫顽场计算代码维护指南

本文面向后续开发者，说明本项目的代码结构、理论背景、数据处理方式、VASP 文件调用关系、公式在代码中的实现位置，以及常见模块修改方法。

主程序包位于 `coercive_field/`。`src/` 是早期探索脚本和辅助脚本，当前命令行工作流不依赖 `src/`。

## 1. 项目目标

本项目从 VASP 的 NEB 路径和 LCALCPOL 极化计算结果中计算外电场下的能量曲线，并估计矫顽场。

核心问题是：

1. 读取每个 NEB image 的零场能量 `U_i`。
2. 读取每个 image 的 LCALCPOL 离子/电子偶极矩。
3. 将 VASP Berry phase 极化放到连续的 polarization branch 上。
4. 用分支修正后的偶极矩变化 `Delta_mu_i` 计算外场修正：

   $$
   \mathbf{E}=E\,\hat{\mathbf{e}},\qquad
   F_i(\mathbf{E})=U_i-0.01\,(\Delta\boldsymbol{\mu}_i\cdot\mathbf{E})
   $$
5. 按 `barrier` 或 `downhill` 判据求使路径判据降到 0 的电场 $E_c$。

单位约定：

- 能量：eV/cell。
- 电场：MV/cm。
- 偶极矩：e*Angstrom。
- 极化：uC/cm^2。
- 体积：Angstrom^3。

## 2. 快速入口

命令行入口有两个：

```text
compute_coercive_field.py
pyproject.toml 中的 console script: compute-coercive-field = coercive_field.cli:main
```

直接运行：

```bash
python compute_coercive_field.py
```

指定路径和外场方向：

```bash
python compute_coercive_field.py --root Pca2_1/xyz-000-011 --field-direction 0,0,-1 --output-dir results_001
```

PowerShell 中负数方向建议用等号，避免参数解析歧义：

```powershell
python compute_coercive_field.py --field-direction=-1,0,-1 --output-dir results_101
```

## 3. 代码模块地图

| 文件                              | 主要作用                                                   | 修改会影响什么                                               |
| --------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------ |
| `compute_coercive_field.py`     | 兼容入口，只调用`coercive_field.cli.main()`              | 通常不用改。若想保持旧命令可用，入口应继续简单。             |
| `coercive_field/cli.py`         | 定义命令行参数、默认样品路径、调用分析和输出               | 改默认路径、默认参数、增加 CLI 参数时改这里。                |
| `coercive_field/models.py`      | 定义`RawImage`、`BranchImage`、`PathResult` 数据结构 | 增加跨模块传递的新数据字段时改这里，并同步读入/输出。        |
| `coercive_field/vasp_io.py`     | 读取 VASP 原始文件：结构、能量、LCALCPOL 偶极矩            | 改文件格式、能量来源、读取规则、路径布局时改这里。           |
| `coercive_field/branch.py`      | Berry branch 选择和连续性修正                              | 改极化分支算法、连续性阈值逻辑、末态分支选择时改这里。       |
| `coercive_field/directions.py`  | 将用户输入方向转为笛卡尔单位矢量                           | 增加方向坐标系或改变晶格方向解释时改这里。                   |
| `coercive_field/constants.py`   | 单位换算常数                                               | 改单位体系或检查换算因子时改这里。                           |
| `coercive_field/calculation.py` | 主计算流程、外场修正、势垒和矫顽场求解                     | 改物理公式、判据、路径范围、矫顽场定义时改这里。             |
| `coercive_field/outputs.py`     | 写 CSV 和画图                                              | 增删输出列、改图像样式、改扫描电场网格时改这里。             |
| `run_all_fe2_paths.sh`          | 批量扫描多个 FE2 路径和多个外场方向                        | 改批量路径匹配、默认方向集合、汇总方式时改这里。             |
| `src/*.py`                      | 早期极化/NEB 辅助脚本                                      | 当前主流程不用。需要复用前先检查脚本底部是否有直接执行代码。 |

## 4. 主调用链

主流程从 CLI 到输出的调用链如下：

```text
compute_coercive_field.py
  -> coercive_field.cli.main()
       -> build_parser()
       -> default_samples() 或 parse_sample_args()
       -> calculation.analyze_path(...) 对每个 sample 计算
            -> vasp_io.read_raw_images(...)
                 -> discover_image_count()
                 -> read_relative_energies()
                      -> parse_oszicar_energy()
                      -> parse_neb_dat()
                 -> image_structure_path()
                 -> parse_structure()
                 -> parse_lcalcpol_dipoles()
                 -> polarization_quanta()
            -> directions.resolve_direction()
            -> branch.build_branch_images(...)
            -> coercive_requirements()
            -> coercive_field_from_requirements()
       -> outputs.write_outputs(...)
            -> 写 CSV
            -> plot_barriers()
            -> plot_profiles()
            -> plot_field_profile_families()
            -> plot_polarization_profiles()
            -> plot_branch_polarization_profiles()
```

数据对象的生命周期：

```text
VASP 文件
  -> RawImage: 原始能量、晶格、体积、偶极矩、原始极化、polarization quantum
  -> BranchImage: Berry 分支编号、修正后极化、修正后偶极矩、连续性诊断
  -> PathResult: 一个 sample 的完整结果、方向、判据、Ec、notes
  -> CSV/PNG 输出
```

## 5. 输入目录和文件依赖

默认根目录是：

```text
Pca2_1/xyz-000-011
```

如果不传 `--sample`，`cli.default_samples()` 会读取 4 个 sample：

```text
undoped-HfO2: <root>
a-HZO:        <root>/a-<root-name>
b-HZO:        <root>/b-<root-name>
c-HZO:        <root>/c-<root-name>
```

例如 `root=Pca2_1/xyz-000-011` 时，默认样品路径是：

```text
Pca2_1/xyz-000-011
Pca2_1/xyz-000-011/a-xyz-000-011
Pca2_1/xyz-000-011/b-xyz-000-011
Pca2_1/xyz-000-011/c-xyz-000-011
```

每个 sample 目录需要包含：

```text
00/OSZICAR
01/OSZICAR
...
08/OSZICAR

00/CONTCAR 或 00/POSCAR
01/CONTCAR 或 01/POSCAR
...
08/CONTCAR 或 08/POSCAR

1-dipol/00/OUTCAR
1-dipol/01/OUTCAR
...
1-dipol/08/OUTCAR
```

可选：

```text
neb.dat
```

## 6. VASP 文件由哪些代码读取

| VASP 文件                           | 读取函数                             | 读取内容                 | 用途                                                  |
| ----------------------------------- | ------------------------------------ | ------------------------ | ----------------------------------------------------- |
| `<sample>/<image>/OSZICAR`        | `vasp_io.parse_oszicar_energy()`   | 最后一个`E0=` 数值     | 构造零场相对能量$U_i=E_i-E_0$。                     |
| `<sample>/neb.dat`                | `vasp_io.parse_neb_dat()`          | 第 3 列相对能量          | `--energy-source neb/auto` 时作为能量来源或校验。   |
| `<sample>/<image>/CONTCAR`        | `vasp_io.parse_structure()`        | scale 和 3 个晶格矢量    | 计算体积`Omega_i`、polarization quantum、方向换算。 |
| `<sample>/<image>/POSCAR`         | `vasp_io.image_structure_path()`   | CONTCAR 缺失时的备用结构 | 与 CONTCAR 用途相同。                                 |
| `<sample>/1-dipol/<image>/OUTCAR` | `vasp_io.parse_lcalcpol_dipoles()` | `p[ion]` 和 `p[elc]` | 得到原始总偶极矩`mu_raw_i`。                        |

注意：

- 能量默认从 NEB image 目录的 `OSZICAR` 读，不从 `1-dipol/` 读。
- 极化只从 `1-dipol/<image>/OUTCAR` 读，因为这些目录应是对 NEB 结构做的 LCALCPOL 静态计算。
- 结构当前从 NEB image 目录读，即 `<sample>/<image>/CONTCAR`，不是从 `1-dipol/<image>/CONTCAR` 读。

## 7. 理论背景和代码公式

### 7.1 NEB 零场相对能量

每个 image 的零场能量记作 $E_i$。代码以 image 00 为参考：

$$
U_i=E_i-E_0
$$

实现位置：

- `vasp_io.read_relative_energies()`
- `RawImage.rel_energy_eV`

`saddle_image` 是 $U_i$ 最大的 image：

$$
i_{\mathrm{saddle}}=\arg\max_i U_i
$$

$$
\Delta U_{\mathrm{barrier}}=U_{i_{\mathrm{saddle}}}
$$

实现位置：

- `calculation.analyze_path()`

### 7.2 LCALCPOL 偶极矩到极化

VASP LCALCPOL 输出离子偶极矩和电子偶极矩：

$$
\boldsymbol{\mu}^{\mathrm{raw}}_i
=\mathbf{p}^{\mathrm{ion}}_i+\mathbf{p}^{\mathrm{elc}}_i
$$

把偶极矩密度换算成极化：

$$
\mathbf{P}^{\mathrm{raw}}_i
=C_P\,\frac{\boldsymbol{\mu}^{\mathrm{raw}}_i}{\Omega_i}
$$

其中：

$$
C_P=\texttt{UC\_CM2\_PER\_EANG\_A3}=1602.176634
$$

实现位置：

- `vasp_io.parse_lcalcpol_dipoles()`
- `vasp_io.read_raw_images()`
- `constants.UC_CM2_PER_EANG_A3`

### 7.3 Polarization quantum

Berry phase 极化只在 polarization quantum 模意义下确定。沿晶格矢量 `a,b,c` 平移一个电子电荷对应：

$$
\mathbf{Q}_a=C_P\,\frac{\mathbf{a}}{\Omega},\qquad
\mathbf{Q}_b=C_P\,\frac{\mathbf{b}}{\Omega},\qquad
\mathbf{Q}_c=C_P\,\frac{\mathbf{c}}{\Omega}
$$

任意等价分支：

$$
\mathbf{P}^{\mathrm{branch}}_i
=\mathbf{P}^{\mathrm{raw}}_i
+n_{i1}\mathbf{Q}_{a,i}
+n_{i2}\mathbf{Q}_{b,i}
+n_{i3}\mathbf{Q}_{c,i}
$$

实现位置：

- `vasp_io.polarization_quanta()`
- `branch.direct_branch_to_reference()`
- `branch.build_branch_images()`

### 7.4 Berry branch 连续性选择

默认 `--branch-mode center-zero`：

1. 从 image 00 开始。
2. 对后续 image 搜索整数分支 `n_i`，使当前 `P_branch_i` 与前一个 image 的 `P_branch_{i-1}` 连续。
3. 对整条已经连续的路径加同一个整数平移，使几何中点 image 的极化最接近 0。
4. 对非 180 度末态，如果可推断末态方向，则用末态参考只修正右半段分支歧义，不强迫整条路径线性。

连续分支选择不再穷举附近整数分支，而是在 polarization quantum 基底中
直接计算需要平移的整数倍。对第 $i$ 个 image，以上一个已经修正好的
$\mathbf{P}^{\mathrm{branch}}_{i-1}$ 为参考：

$$
Q_i^T\mathbf{c}_i
=\mathbf{P}^{\mathrm{branch}}_{i-1}-\mathbf{P}^{\mathrm{raw}}_i
$$

$$
\mathbf{n}_i=\lfloor \mathbf{c}_i+0.5\rfloor
$$

$$
\mathbf{P}^{\mathrm{branch}}_i
=\mathbf{P}^{\mathrm{raw}}_i+\mathbf{n}_i Q_i
$$

这里 $Q_i$ 的三行是 $\mathbf{Q}_{a,i},\mathbf{Q}_{b,i},\mathbf{Q}_{c,i}$。
对于非正交晶胞，不能把 Cartesian 分量分别除以某个标量 quantum；
必须先解上面的 quantum-basis 系数。使用 $\lfloor c+0.5\rfloor$ 是为了取
最接近参考极化的整数分支；如果目标只是把单个极化值折回 $[0,Q)$ 主区间，
才会使用简单的 $\lfloor c\rfloor$。

实现位置：

- `branch.fold_continuous_from_previous()`
- `branch.fold_to_target()`
- `branch.refine_forward_branch_to_endpoint()`
- `branch.build_branch_images()`

输出诊断：

- `branch_correction.csv`
- `polarization_branches.csv`
- `branch_polarization_profiles.png`

### 7.5 外场方向和分支诊断方向

用户输入的 `field_direction` 默认按晶轴 $(a,b,c)$ 坐标解释。例如
`-1,0,0` 表示沿 $-\mathbf{a}$ 方向施加电场。代码先用第一张结构的完整
晶格矢量把这个晶轴方向转换成笛卡尔单位矢量 $\hat{\mathbf{e}}$，再用于
能量修正：

$$
\Delta\mu_i^{\parallel}
=\Delta\boldsymbol{\mu}_i\cdot\hat{\mathbf{e}}
$$

分支方向 `branch_direction` 主要用于分支投影和单调性诊断，不等价于外场方向。

方向输入支持 3 种坐标解释：

$$
\text{cartesian:}\quad \mathbf{r}=\mathbf{d}
$$

$$
\text{lattice:}\quad \mathbf{r}=d_1\mathbf{a}+d_2\mathbf{b}+d_3\mathbf{c}
$$

$$
\text{lattice-diagonal:}\quad
\mathbf{r}=(d_1a_x,\ d_2b_y,\ d_3c_z)
$$

$$
\hat{\mathbf{e}}=\frac{\mathbf{r}}{|\mathbf{r}|}
$$

实现位置：

- `directions.parse_vector()`
- `directions.resolve_direction()`

默认 `--field-coordinates lattice`，即按完整晶格矢量解释
`--field-direction`。`lattice-diagonal` 仅保留为旧 Excel 兼容模式，只按
$a_x,b_y,c_z$ 三个对角长度近似归一化；非正交晶胞不应默认使用。

### 7.6 分支修正偶极矩变化

外场修正不用线性目标极化，而是用 VASP LCALCPOL 数据经 Berry branch 修正后的偶极矩：

$$
\boldsymbol{\mu}^{\mathrm{branch}}_i
=\mathbf{P}^{\mathrm{branch}}_i\,\frac{\Omega_i}{C_P}
$$

$$
\Delta\boldsymbol{\mu}_i
=\boldsymbol{\mu}^{\mathrm{branch}}_i
-\boldsymbol{\mu}^{\mathrm{branch}}_0
$$

实现位置：

- `branch.build_branch_images()`
- `calculation.branch_delta_dipoles()`

### 7.7 外场能量修正

电场下相对自由能应写成矢量点乘。令

$$
\mathbf{E}=E\,\hat{\mathbf{e}}
$$

则

$$
F_i(\mathbf{E})=U_i-C_E\,(\Delta\boldsymbol{\mu}_i\cdot\mathbf{E})
$$

其中：

$$
C_E=\texttt{EV\_PER\_EANG\_MV\_CM}=0.01
$$

物理含义：

$$
1\,(e\cdot\mathrm{\AA})\times1\,(\mathrm{MV/cm})=0.01\,\mathrm{eV}
$$

实现位置：

- `constants.EV_PER_EANG_MV_CM`
- `calculation.corrected_relative_energies()`

输出位置：

- `corrected_profiles.csv`
- `energy_profiles.png`
- `field_profiles.png`

### 7.8 矫顽场判据

代码支持两种判据。

`barrier` 判据：

$$
B(E)=\max_i F_i(E)-F_0(E)
$$

$$
B(E_c)=0
$$

`downhill` 判据：

$$
D(E)=\max_i\left[F_i(E)-F_{i-1}(E)\right]
$$

$$
D(E_c)=0
$$

路径 image 范围由 `--barrier-mode` 决定：

```text
to-saddle:     i = 1..saddle_image
intermediate:  i = 1..N-2
full-path:     i = 1..N-1
```

实现位置：

- `calculation.barrier_indices()`
- `calculation.path_value_at_field()`
- `calculation.candidate_indices_for_mode()`
- `calculation.coercive_requirements()`

### 7.9 解析式求 Ec

每个约束都可以写成：

$$
\Delta E-0.01\,E\,\Delta\mu^{\parallel}\le 0
$$

若 $\Delta E>0$ 且 $\Delta\mu^{\parallel}>0$：

$$
E_{\mathrm{required}}
=\frac{\Delta E}{0.01\,\Delta\mu^{\parallel}}
$$

最终：

$$
E_c=\max_i E_{\mathrm{required},i}
$$

如果某个正能量约束 `Delta_mu_parallel <= 0`，说明当前外场方向无法降低该约束，代码返回 `inf` 并在 `notes` 中记录。

实现位置：

- `calculation.coercive_requirements()`
- `calculation.coercive_field_from_requirements()`

输出位置：

- `coercive_field_summary.csv`
- `field_projection.csv`
- `barrier_vs_field.csv`
- `barrier_vs_field.png`

## 8. 数据处理流程细节

### 8.1 读取原始数据

`vasp_io.read_raw_images(label, path, energy_source)` 返回：

```text
list[RawImage], notes
```

每个 `RawImage` 含有：

- sample 标签和路径。
- image 编号。
- 能量文件路径、结构文件路径、dipole OUTCAR 路径。
- `energy_eV` 和 `rel_energy_eV`。
- 晶格矩阵 `lattice_A` 和体积 `volume_A3`。
- LCALCPOL 的 `ion_dipole_eA`、`electronic_dipole_eA`。
- `raw_dipole_eA = ion + electronic`。
- `raw_polarization_uC_cm2`。
- `quantum_vectors_uC_cm2`。

### 8.2 修正 Berry branch

`branch.build_branch_images(...)` 接收 `raw_images`，返回 `list[BranchImage]`。

每个 `BranchImage` 含有：

- `branch_index = [n1, n2, n3]`。
- `branch_polarization_uC_cm2`。
- `branch_dipole_eA`。
- 相邻 image 的极化变化。
- `max_abs_component_jump_from_prev_uC_cm2`。
- `component_continuity_ok`。
- `monotonic_along_branch`。

这一步是整个计算中最容易影响结果的部分。后续外场修正完全依赖 `branch_dipole_eA`。

### 8.3 求解单个 sample 的结果

`calculation.analyze_path(...)` 把一个 sample 从原始文件处理到 `PathResult`。

它做的事情：

1. 调 `read_raw_images()` 读 VASP 文件。
2. 用第一张结构的晶格解析外场方向和 branch 方向。
3. 推断或解析末态分支参考方向。
4. 调 `build_branch_images()` 生成连续分支。
5. 计算 `projected_delta`，即 $(\boldsymbol{\mu}^{\mathrm{branch}}_i-\boldsymbol{\mu}^{\mathrm{branch}}_0)\cdot\hat{\mathbf{e}}$。这是解析求 $E_c$ 时未知标量场强 $E$ 的线性系数；实际能量修正使用 $\Delta\boldsymbol{\mu}_i\cdot\mathbf{E}$。
6. 按判据求 `requirements`。
7. 得到 `coercive_field_MV_cm` 和 `limiting_image`。
8. 把结果打包成 `PathResult`。

### 8.4 写输出

`outputs.write_outputs(results, output_dir, field_step, field_max)` 接收多个 sample 的 `PathResult`，生成 CSV 和图片。

输出文件：

```text
coercive_field_summary.csv
raw_vasp_data.csv
branch_correction.csv
polarization_branches.csv
field_projection.csv
barrier_vs_field.csv
corrected_profiles.csv
polarization_profiles.csv
barrier_vs_field.png
energy_profiles.png
field_profiles.png
polarization_profiles.png
branch_polarization_profiles.png
```

`field_step` 和 `field_max` 只影响 `barrier_vs_field.csv/png` 的扫描网格，不影响解析求得的 `E_c`。

## 9. 输出表和上游来源

| 输出文件                       | 主要来源                                        | 作用                                                 |
| ------------------------------ | ----------------------------------------------- | ---------------------------------------------------- |
| `coercive_field_summary.csv` | `PathResult`                                  | 每个 sample 一行，最终结果总表。                     |
| `raw_vasp_data.csv`          | `RawImage`                                    | 检查 VASP 原始读取和单位换算。                       |
| `branch_correction.csv`      | `BranchImage`                                 | 检查 Berry branch 整数平移和连续性。                 |
| `polarization_branches.csv`  | `BranchImage`                                 | 与`branch_correction.csv` 相同，保留兼容旧文件名。 |
| `field_projection.csv`       | `RawImage` + `BranchImage` + `PathResult` | 检查外场投影、每个约束的 required field。            |
| `barrier_vs_field.csv`       | `barrier_at_field()`                          | 电场扫描曲线，画`barrier_vs_field.png`。           |
| `corrected_profiles.csv`     | `corrected_relative_energies()`               | 只写 0 和 Ec 下的能量曲线。                          |
| `polarization_profiles.csv`  | `target_polarization()` + branch 极化         | 检查目标路径和实际分支路径。                         |

## 10. 常见修改指南

### 10.1 增加一个新的命令行参数

修改：

1. `coercive_field/cli.py` 的 `build_parser()` 添加参数。
2. `cli.main()` 调用 `analyze_path()` 或 `write_outputs()` 时传入。
3. 如果参数需要长期保存在结果中，给 `PathResult` 增加字段。
4. 如果需要输出，修改 `outputs.write_outputs()` 的对应 CSV 行。

检查：

```bash
python compute_coercive_field.py --help
python compute_coercive_field.py <新参数> --output-dir results_test
```

### 10.2 改默认 sample 发现逻辑

修改：

- `coercive_field/cli.py` 中的 `default_samples(root)`。

当前逻辑假设掺杂目录在 root 内部：

```text
root/a-<root-name>
root/b-<root-name>
root/c-<root-name>
```

如果以后目录变成同级目录或更多掺杂类型，只改这个函数通常就够。

### 10.3 改能量读取来源

如果只是选择 OSZICAR/neb.dat：

- 用 CLI 参数 `--energy-source oszicar|neb|auto`。

如果要新增 OUTCAR 能量、data-barrier 或其他格式：

1. 在 `vasp_io.py` 添加解析函数。
2. 修改 `read_relative_energies()` 支持新的 `source`。
3. 在 `cli.py` 的 `--energy-source choices` 增加新选项。
4. 保证返回的是相对 image 00 的 `rel_energies`。

注意：后续计算假设 `rel_energies[0] == 0`。

### 10.4 改结构文件来源

当前结构优先：

```text
<sample>/<image>/CONTCAR
<sample>/<image>/POSCAR
```

修改位置：

- `vasp_io.image_structure_path()`

如果想改成读 `1-dipol/<image>/CONTCAR`，应同时确认 LCALCPOL 静态计算是否真的保留了同一结构和晶格。体积 `Omega_i` 会影响极化换算、polarization quantum 和 branch dipole。

### 10.5 改 LCALCPOL OUTCAR 路径或格式

当前读取：

```text
<sample>/1-dipol/<image>/OUTCAR
```

修改位置：

- `vasp_io.read_raw_images()` 中构造 `path / "1-dipol" / f"{image:02d}" / "OUTCAR"` 的地方。
- 如果 VASP 输出文本格式变了，改 `vasp_io.parse_lcalcpol_dipoles()` 的正则表达式。

### 10.6 改 Berry branch 算法

优先判断你要改哪一层：

| 目标                 | 修改位置                                                                            |
| -------------------- | ----------------------------------------------------------------------------------- |
| 改默认分支模式       | `cli.py` 的 `--branch-mode` 默认值。                                            |
| 增加新的 branch mode | `branch.build_branch_images()` 加分支，并在 `cli.py` choices 加选项。           |
| 改连续分支取整规则   | `branch.direct_branch_to_reference()` 和 `branch.fold_continuous_from_previous()`。 |
| 改目标路径弱约束     | `branch.fold_continuous_near_target()`。                                          |
| 改中点置零逻辑       | `branch.build_branch_images()` 中 `global_shift` 部分。                         |
| 改非 180 度末态修正  | `branch.endpoint_reference()` 和 `branch.refine_forward_branch_to_endpoint()`。 |
| 改连续性警告阈值     | CLI 参数`--continuity-threshold` 或 `build_branch_images()` 里的判断。          |

改完后重点检查：

- `branch_correction.csv` 的 `branch_n1/2/3`。
- `max_abs_component_jump_from_prev_uC_cm2`。
- `component_continuity_ok`。
- `branch_polarization_profiles.png` 是否出现不合理跳变。

### 10.7 改外场方向解释

修改：

- `directions.resolve_direction()`
- `cli.py` 的 `--field-coordinates` 或 `--branch-coordinates` choices。

注意：

- `field_direction` 会进入能量公式。
- `branch_direction` 主要是诊断和少量 branch 参考逻辑。
- 不要把二者混在一起，除非明确希望 branch 诊断方向等于外场方向。

### 10.8 改外场能量修正公式

当前公式应按矢量点乘理解：

$$
\mathbf{E}=E\,\hat{\mathbf{e}},\qquad
F_i(\mathbf{E})=U_i-0.01\,(\Delta\boldsymbol{\mu}_i\cdot\mathbf{E})
$$

修改位置：

- `calculation.corrected_relative_energies()`
- `calculation.coercive_requirements()`
- `outputs.write_outputs()` 中 `slope`、`constraint_slope` 和 `required_field` 的计算。

如果只改常数：

- `constants.EV_PER_EANG_MV_CM`

注意：如果公式改了但 `coercive_requirements()` 仍用旧解析式，最终 `E_c` 会和扫描曲线不一致。因此能量修正和解析求解必须一起改。

### 10.9 增加新的矫顽场判据

修改步骤：

1. 在 `cli.py` 的 `--coercive-criterion choices` 增加新名字。
2. 在 `calculation.path_value_at_field()` 增加该判据的路径值。
3. 在 `calculation.coercive_requirements()` 增加解析求解逻辑。
4. 在 `outputs.write_outputs()` 中增加每行约束解释和 `required_field` 输出。
5. 修改图的 ylabel：`outputs.plot_barriers()`。

如果新判据不能解析求解，可以：

- 保留 `path_value_at_field()`。
- 新增数值求根函数。
- 在 `analyze_path()` 中根据 criterion 选择解析或数值求解。

### 10.10 改参与计算的 image 范围

已有模式：

```text
to-saddle
intermediate
full-path
```

修改位置：

- `calculation.barrier_indices()`
- `calculation.candidate_indices_for_mode()`
- `cli.py` 的 `--barrier-mode choices`

这两个函数要保持一致：

- `barrier_indices()` 给扫描/画图的路径值用。
- `candidate_indices_for_mode()` 给解析求 `E_c` 用。

### 10.11 增删输出列

修改：

- `outputs.write_outputs()`

常用辅助函数：

- `vector_columns(prefix, vec)` 展开三维矢量为 `prefix_x/y/z`。
- `matrix_columns(prefix, mat)` 展开 3x3 矩阵为 `prefix_a_x` 等列。

如果新增字段来自新数据结构：

1. 先在 `models.py` 添加字段。
2. 在生成 `RawImage`、`BranchImage` 或 `PathResult` 的地方填入。
3. 在 `outputs.py` 写入 CSV。

### 10.12 改图像样式

修改：

- `outputs.PLOT_STYLE`
- `outputs.SERIES_COLORS`
- `outputs.SERIES_MARKERS`
- 各 `plot_*()` 函数。

当前 Matplotlib 设置启用了：

```python
"text.usetex": True
```

如果运行环境没有 LaTeX，画图可能失败。可以在 `PLOT_STYLE` 中改为：

```python
"text.usetex": False
```

### 10.13 改批量扫描方向

修改：

- `run_all_fe2_paths.sh` 中 `FIELD_DIRECTIONS` 数组。

当前包含：

```text
001:  0,0,-1
010:  0,1,0
101: -1,0,-1
m110:-1,1,0
111: -1,-1,-1
```

批量脚本会扫描 phase 根目录下一层中名字匹配：

```text
[xyz][xyz][xyz]-[01][01][01]-[01][01][01]
```

的路径目录。

### 10.14 使用 `src/` 中的旧脚本

`src/` 主要是早期探索代码：

| 文件                                      | 内容                                               |
| ----------------------------------------- | -------------------------------------------------- |
| `src/calculate_polarization.py`         | 单个目录的 LCALCPOL 极化读取和 branch shift 计算。 |
| `src/calculate_polarization_quantum.py` | 从 CONTCAR/POSCAR 计算 polarization quantum。      |
| `src/polarization_neb.py`               | 旧版 NEB 极化、能量、外场扫描流程。                |
| `src/utils.py`                          | 当前基本为空。                                     |

注意：`src/polarization_neb.py` 文件底部有直接执行示例代码，导入该文件时可能立即运行扫描。因此把里面函数复用到主包前，建议先把示例代码放进：

```python
if __name__ == "__main__":
    ...
```

当前主包 `coercive_field/` 已经重新实现了这些功能，后续开发优先改主包。

## 11. 数据一致性检查

修改后建议至少做以下检查：

```bash
python compute_coercive_field.py --root Pca2_1/xyz-000-011 --output-dir results_dev_check
```

然后检查：

1. `results_dev_check/coercive_field_summary.csv`

   - `coercive_field_MV_cm` 是否有限或 notes 是否合理。
   - `limiting_constraint` 是否符合预期。
2. `results_dev_check/branch_correction.csv`

   - `component_continuity_ok` 是否大多为 `True`。
   - `max_abs_component_jump_from_prev_uC_cm2` 是否异常大。
3. `results_dev_check/field_projection.csv`

   - 限制约束的 `required_field_MV_cm` 是否等于 summary 中的 `E_c`。
   - 如果 `can_be_lowered_by_field=False`，确认外场方向是否选反。
4. 图片：

   - `branch_polarization_profiles.png` 看 branch 是否跳跃。
   - `barrier_vs_field.png` 看判据是否随场下降到 0。
   - `energy_profiles.png` 看 Ec 处是否达到所选判据。

## 12. 新开发的建议顺序

当需求不确定时，建议按下面顺序修改：

1. 先确认要改的是输入、branch、外场公式、判据还是输出。
2. 找到上面对应的最小模块。
3. 如果新增跨模块数据，先改 `models.py`。
4. 再改产生数据的模块。
5. 再改 `outputs.py` 写出诊断列。
6. 最后改 CLI 默认值或参数。
7. 用一个小输出目录跑一遍并检查 CSV 中间量。

原则是：所有影响最终 `E_c` 的新中间量，都应该能在 CSV 中追踪。这样后续排查时可以从 `summary -> field_projection -> branch_correction -> raw_vasp_data` 逐层回溯。
