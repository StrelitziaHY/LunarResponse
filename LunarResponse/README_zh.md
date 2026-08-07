# LunarResponse（月球响应函数计算）

[English README](README.md)

**从内部结构模型到月表响应曲线，一条命令完成。**
输入月球（或行星）内部结构模型，本程序自动完成：模型插值加密 → 用修改版
MINEOS 计算自由振荡本征模式 → 模式叠加得到随频率变化的月表位移响应函数，
包括 A 类（Dyson 体力驱动）与 B 类（引潮位驱动）响应，适用于引力波 /
引潮响应研究。

```bash
python lunar_response.py all --model examples/MoonModelSimp.txt --name demo --out run_demo
```

这一条命令即可在 `run_demo/` 中生成模式表与响应曲线
（`response_functions.txt` + 两张 PNG 图）。

---

## 功能特点

- **端到端流程**：简化分层模型 → 加密 MINEOS 卡片模型 → 本征频率与本征函数
  （球型与环型）→ 响应曲线。各阶段也可单独运行（`model` / `eigen` / `response`）。
- **扩容版 MINEOS**：径向层数上限从原版 `mk=350` 提高到 `mk=50000`，
  可求解插值到数万层的模型（已用 `mk=60000` 验证 55001 层）。物理算法未改动，
  结果与原版 MINEOS 逐位一致。
- **预编译 Windows 可执行文件**（`bin/`）：日常使用无需安装编译器。
- **稳健的模式表解析**（不再依赖"删除固定行数"），数值计算全向量化
  （典型 l=2 算例约 350 个模式，全流程约 30 秒）。
- **流体层**（如液态核）自动识别，也可显式指定。

## 运行环境

- Python 3 + `numpy`（绘图需要 `matplotlib`）。
- Windows：直接使用 `bin/` 中的预编译程序即可。
- 其他平台、或需要修改层数上限：需要 `gfortran`（Windows 下可运行
  `python install_mingw.py` 安装便携版编译器，无需管理员权限）。

## 快速上手

```bash
# 全流程（默认参数即复现均匀月球 l=2 基准算例）
python lunar_response.py all --model examples/MoonModelSimp.txt --name demo --out run_demo

# 分阶段运行
python lunar_response.py model    --model examples/MoonModelSimp.txt --name demo --out run_demo
python lunar_response.py eigen    --model-file run_demo/demo.txt --out run_demo
python lunar_response.py response --folder run_demo --code demo_code.txt
```

程序自动定位 MINEOS 可执行文件：`--exe-dir` 参数 > 环境变量 `LUNAR_MINEOS_DIR`
> 脚本同级 `bin/` 目录。

### 输入模型格式

六列，从**表面到球心**排列：

```
r(km)  vp(km/s)  vs(km/s)  rho(g/cm^3)  Q_kappa  Q_mu
```

流体层令 `vs = 0` 即可（程序自动识别流体块并写入 MINEOS 卡头的 `nic/noc`）。

## 主要可调参数

| 参数 | 含义 | 默认值 |
|---|---|---|
| `--dr-km` / `--layers` | 插值精度：目标最大层厚（km）**或**目标总层数（二选一） | legacy 方案 `--mult/--mj/--crit` |
| `--nic`, `--noc` | 流体层界（内核固体侧 / 核幔边界流体侧层号） | 自动识别（`vs=0, rho>0`） |
| `--nmin`, `--nmax` | 径向阶数范围（可不从最小 n 开始） | 0, 400 |
| `--wmin`, `--wmax` | 频率范围（mHz） | 0.1, 300 |
| `--wgrav` | 自引力阈值（mHz），高于此频率忽略重力项（约快 3 倍） | 0.5 |
| `--jcom` | 模式类型：1 径向，2 环型，3 球型，4 内核环型 | 3 |
| `--lmin`, `--lmax` | 角阶数范围 | 2, 2 |
| `--eps` | Runge–Kutta 积分精度（环型模式不宜小于 ~2e-4） | 1e-8 |
| `--save-format` | `compact`（5 列）或 `full`（7 列，含未减静态的 B 类） | compact |
| `--static-rad`, `--static-hor` | B 类静态极限系数（×R；l=2 时为 0.5 / 0.25） | 0.5, 0.25 |

完整参数列表见 `python lunar_response.py <子命令> -h`。

## 输出约定：A 类与 B 类响应

- **A 类 `T_A`**：直接体力（Dyson forcing）驱动下的表面位移响应；
- **B 类 `T_B`**：引潮位驱动下的表面位移响应。其静态极限为 `0.5·R`（径向）、
  `0.25·R`（切向，l=2）；动力学部分 `|T_B − 静态极限|` 与 A 类量级相近。

`compact` 格式 5 列：`freq, T_A_r, T_A_h, |T_B_r−0.5R|, |T_B_h−0.25R|`；
`full` 格式额外包含未减静态的 `|T_B_r|`、`|T_B_h|`。输出文件头中写有明确定义。

## 重新编译 MINEOS / 调整层数上限

```bash
python install_mingw.py             # 一次性：无 gfortran 时安装便携编译器（Windows）
python build_mineos.py              # 默认 mk=50000 编译到 ./bin
python build_mineos.py --mk 100000  # 一条命令提高层数上限
```

`build_mineos.py` 在临时副本中统一替换全部 `parameter (mk=...)`（`src/`
源码永不修改），用 `-O2` 半静态链接编译，自动附带运行时 DLL 并做冒烟测试。
**注意**：模型层数超过 mk 时程序会**挂起而非报错**，务必保证 层数 < mk。

## 设计基础

- **模式叠加法**。本征频率与本征函数由 MINEOS 计算（`minos_bran` 打靶法 +
  Runge–Kutta 积分；`eigcon`；`eigen2asc`），归一化约定见 Woodhouse & Dahlen
  (1978)。模式耦合系数（引潮载荷与 Dyson 载荷，后者对剪切模量梯度做分部积分）
  在模型网格上用 Simpson 积分求值，再以阻尼振子传播子
  `1/(ω_n² − ω² + iωω_n/Q_n)` 对所有模式求和。
- **相对上游 MINEOS 的源码修改**（刻意保持最小）：三个文件中
  `mk=350 → 50000`；`fdb_eigen` 辅助子程序内联进 `eigcon.f` 以单文件编译；
  注释掉一处 `loctime()` 调用。数值结果与上游完全一致。
- 已对原 notebook 流程做逐位验证（356 个球型模式、567 个环型模式及最终
  响应曲线全部精确一致）。

## 目录结构

```
lunar_response.py    端到端主程序（模型 → 本征模式 → 响应曲线）
build_mineos.py      一键构建 MINEOS（层数上限 mk 可调）
install_mingw.py     便携 gfortran 安装器（Windows，免管理员权限）
src/                 修改版 MINEOS Fortran 源码（见"设计基础"）
bin/                 预编译 Windows 可执行文件 + 运行时 DLL
examples/            示例输入模型（均匀月球）
LICENSE              GPL v2（继承自 MINEOS）
```

## 许可证

Fortran 源码派生自 MINEOS，遵循 **GNU 通用公共许可证 v2**（见 `LICENSE`）；
Python 驱动程序按相同条款分发。
