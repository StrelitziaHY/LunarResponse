# -*- coding: utf-8 -*-
"""
lunar_response.py — 月球响应函数端到端计算程序
=================================================

用一个输入模型文件，一条龙完成：

    简化模型 (r_km, vp, vs, rho, Qk, Qs)
        -> 1. 插值加密并生成 mineos 卡片式模型
        -> 2. 调用 mineos 三步程序 (1calc / 2caltoeigen / 3eigentoasc) 计算本征模式
        -> 3. 稳健解析模式表与本征函数，计算模式积分
        -> 4. 叠加模式得到响应函数曲线 T_A / T_B，保存数据与图

对比原 notebook 流程 (EigCalc_UnityVer / TransfCalc_2025new) 的改进：
  * 单一脚本 + 命令行，输入初始模型直接输出响应曲线；
  * 插值、积分、响应叠加全部向量化（原代码为逐点 Python 循环）；
  * res1 模式表不再靠"删除前 Nc+9/Nc+11 行"这种脆弱定位，而是按表头自动定位、
    自动剥离 's'/'t' 后缀与 '***' 标记；
  * 不移动/删除任何输入文件，所有产物写入独立输出目录；
  * 球型 (jcom=3) 与环型 (jcom=2) 模式统一处理。

用法示例：
  # 全流程（需要三个 mineos 可执行文件在当前目录或用 --exe-dir 指定）
  python lunar_response.py all --model MoonModelSimp.txt --name run1 \
      --mult 8000 --mj 200 --jcom 3 --lmin 2 --lmax 2 --nmax 400 --wmax 300

  # 只做模型加密
  python lunar_response.py model --model MoonModelSimp.txt --name RM25_test

  # 已有模式库，只算响应
  python lunar_response.py response --folder RM25_8001 --code RM25_8001_code.txt
"""

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- 工具


def simpson_irregular(y, x):
    """非均匀网格 Simpson 积分（与 scipy.integrate.simps 等价的三点抛物线公式）。

    样本数为偶数时，最后一个区间用梯形补充。
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 3:
        return float(np.trapezoid(y, x))
    result = 0.0
    # 每次跨两点（i, i+1, i+2）
    last = n - 1 if (n % 2 == 1) else n - 2
    for i in range(0, last - 1, 2):
        h0 = x[i + 1] - x[i]
        h1 = x[i + 2] - x[i + 1]
        hsum = h0 + h1
        hprod = h0 * h1
        hdiv = h0 / h1
        result += (hsum / 6.0) * (
            y[i] * (2.0 - 1.0 / hdiv)
            + y[i + 1] * hsum * hsum / hprod
            + y[i + 2] * (2.0 - hdiv)
        )
    if n % 2 == 0:  # 偶数样本，最后一个区间梯形
        result += 0.5 * (x[-1] - x[-2]) * (y[-1] + y[-2])
    return float(result)


# ---------------------------------------------------------------- 配置


@dataclass
class EigenConfig:
    """mineos 本征模式计算参数（对应 1calc_levels 的交互输入）。"""
    jcom: int = 3          # 1 径向, 2 环型, 3 球型, 4 内核环型
    eps: float = 1e-8      # Runge-Kutta 积分精度（球型建议 1e-8；环型不宜太小，~2e-4）
    wgrav: float = 0.5     # mHz，高于此频率忽略重力项加速计算
    lmin: int = 2
    lmax: int = 2
    wmin: float = 0.1      # mHz
    wmax: float = 300.0    # mHz
    nmin: int = 0
    nmax: int = 400
    dmax_km: float = 6371.0  # 2caltoeigen 的 max depth（>天体半径即取全域）


@dataclass
class ModelConfig:
    """模型插值加密参数。

    三种方案（scheme）：
      "dr"     — 单一精度参数 dr_km：每层细分到层厚 <= dr_km（推荐，精度含义直观）
      "layers" — 单一精度参数 target_layers：按层厚比例分配，总格点数 = target_layers
      "legacy" — 兼容旧 notebook：厚层(>crit_dis_km) mult 份、薄层 mj 份
    """
    scheme: str = "dr"
    dr_km: float = 0.22
    target_layers: int = 8001
    mult: int = 8000
    mj: int = 200
    crit_dis_km: float = 30.0


# ---------------------------------------------------------------- 1. 模型生成


def load_simple_model(path):
    """读取简化模型文件：列 = r(km), vp(km/s), vs(km/s), rho(g/cm3), Qk, Qs，
    从表面到中心排列。返回中心->表面、SI 单位的数组。"""
    m = np.loadtxt(path)
    if m.ndim == 1:
        m = m[None, :]
    m = m[::-1]  # 反转为 中心->表面
    return dict(
        r=m[:, 0] * 1000.0,
        vp=m[:, 1] * 1000.0,
        vs=m[:, 2] * 1000.0,
        rho=m[:, 3] * 1000.0,
        qk=m[:, 4].astype(float),
        qs=m[:, 5].astype(float),
    )


def refine_model(simple, cfg: ModelConfig):
    """按选定方案加密并线性插值（向量化实现；legacy 方案等价于原 notebook 的
    round(a+(b-a)*j/tot, 2) 逐点写法）。返回 9 列 mineos 数组：
    r(m), rho, vp, vs, Qk, Qs, vph, vsh, eta。"""
    r = simple["r"]
    props = [simple[k] for k in ("rho", "vp", "vs", "qk", "qs")]
    dr = np.abs(np.diff(r))

    if cfg.scheme == "dr":
        counts = np.maximum(1, np.ceil(dr / (cfg.dr_km * 1000.0)).astype(int))
    elif cfg.scheme == "layers":
        counts = np.maximum(
            1, np.rint((cfg.target_layers - 1) * dr / dr.sum()).astype(int))
    elif cfg.scheme == "legacy":
        counts = np.where(dr > cfg.crit_dis_km * 1000.0, cfg.mult, cfg.mj)
    else:
        raise ValueError(f"未知插值方案: {cfg.scheme}")

    # 每层生成 counts[i] 个内部插值点（j/tot, j=0..counts-1）
    frac = np.concatenate([np.arange(c) / c for c in counts])
    seg = np.repeat(np.arange(len(counts)), counts)
    t = frac

    r_new = r[seg] + (r[seg + 1] - r[seg]) * t
    cols = [np.round(r_new, 2)]
    for p in props:
        cols.append(np.round(p[seg] + (p[seg + 1] - p[seg]) * t, 2))

    # 追加表面点
    surf = [round(float(r[-1]), 2)] + [round(float(p[-1]), 2) for p in props]
    data = np.vstack([np.column_stack(cols), np.array(surf)])

    # mineos 9 列: r, rho, vp, vs, qk, qs, vph, vsh, eta
    out = np.column_stack([
        data[:, 0], data[:, 1], data[:, 2], data[:, 3],
        data[:, 4], data[:, 5], data[:, 2], data[:, 3],
        np.ones(len(data)),
    ])
    return out


def detect_fluid_indices(arr9):
    """从 9 列加密模型自动识别流体层，返回 (nic, noc)。

    流体层判据：vs == 0 且 rho > 0。取最深（最靠球心）的连续流体块：
      nic = 该块下方固体层数（= 最内固体边界下方层数，无固体内核则为 0）
      noc = 该块最外层流体层的 1-based 层号
    无流体层返回 (0, 0)；若顶部还有海水层等第二流体块，打印警告。
    """
    vs, rho = arr9[:, 3], arr9[:, 1]
    fluid = (vs == 0.0) & (rho > 0.0)
    if not fluid.any():
        return 0, 0
    # 找连续流体块
    idx = np.flatnonzero(fluid)
    breaks = np.where(np.diff(idx) > 1)[0]
    blocks = np.split(idx, breaks + 1)
    deepest = blocks[0]  # 格点从球心到表面，首个块最深
    if len(blocks) > 1:
        print(f"警告: 检测到 {len(blocks)} 个流体块（可能含顶部海水层），"
              f"nic/noc 仅按最深的核部流体块设置。")
    s, e = int(deepest[0]), int(deepest[-1])
    return s, e + 1  # nic=s（下方固体层数），noc=e+1（1-based 最外流体层号）


def write_mineos_models(arr9, prefix, out_dir, nic=None, noc=None):
    """写出 <prefix>_code.txt（纯数据）与 <prefix>.txt（mineos 卡片式模型）。
    nic/noc 为 None 时自动检测流体层。"""
    if nic is None or noc is None:
        nic_auto, noc_auto = detect_fluid_indices(arr9)
        nic = nic_auto if nic is None else nic
        noc = noc_auto if noc is None else noc
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    code_path = out_dir / f"{prefix}_code.txt"
    fmt = lambda row: " ".join(str(v) for v in row)
    with open(code_path, "w") as f:
        f.write("\n".join(fmt(row) for row in arr9) + "\n")

    model_path = out_dir / f"{prefix}.txt"
    with open(model_path, "w") as f:
        f.write(f"{prefix} Model: (  0.0000   0.000) isotropic case \n")
        f.write("  0  1  1 \n")
        f.write(f"  {len(arr9)}  {nic}  {noc} \n")
        f.write("\n".join(fmt(row) for row in arr9) + "\n")
    return model_path, code_path


# ---------------------------------------------------------------- 2. mineos 本征模式计算


def _feed(exe, lines, cwd):
    """向交互式 Fortran 程序按行喂 stdin。"""
    proc = subprocess.Popen(
        [str(exe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, cwd=cwd,
    )
    out, _ = proc.communicate("\n".join(lines) + "\n")
    if proc.returncode != 0:
        raise RuntimeError(f"{exe} 返回码 {proc.returncode}\n{out[-2000:]}")
    return out


def resolve_exe_dir(exe_dir=None):
    """定位 MINEOS 可执行文件目录。搜索顺序：
    1) --exe-dir 参数；2) 环境变量 LUNAR_MINEOS_DIR；
    3) 本脚本同级 mineos_build/bin（build_mineos.py 的默认产物）；
    4) 当前目录。目录中须同时存在三个 exe。
    """
    EXES = ("1calc_levels.exe", "2caltoeigen_levels.exe", "3eigentoasc_levels.exe")
    cands = []
    if exe_dir:
        cands.append(Path(exe_dir))
    env = os.environ.get("LUNAR_MINEOS_DIR")
    if env:
        cands.append(Path(env))
    cands.append(Path(__file__).parent / "bin")
    cands.append(Path(__file__).parent / "mineos_build" / "bin")
    cands.append(Path("."))
    for c in cands:
        if all((c / e).exists() for e in EXES):
            return c.resolve()
    raise FileNotFoundError(
        "未找到 MINEOS 可执行文件（1/2/3*_levels.exe）。"
        "请用 --exe-dir 指定，或先运行 mineos_build/build_mineos.py 编译。")


def run_eigen_pipeline(model_path, cfg: EigenConfig, exe_dir, work_dir,
                       dbname="db"):
    """运行 mineos 三步流程。返回 (mode_table_path, asc_dir)。

    work_dir 下生成 res1 / eig1 / <dbname>.eigen* / <dbname>/ (ASC 文件)。
    所有中间文件保留在 work_dir，不挪动输入模型。
    """
    exe_dir = resolve_exe_dir(exe_dir)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    model_path = Path(model_path).resolve()

    # --- 1. minos_bran：算本征频率 ---
    _feed(
        exe_dir / "1calc_levels.exe",
        [str(model_path), "res1", "eig1",
         f"{cfg.eps} {cfg.wgrav}", str(cfg.jcom),
         f"{cfg.lmin} {cfg.lmax} {cfg.wmin} {cfg.wmax} {cfg.nmin} {cfg.nmax}"],
        cwd=work_dir,
    )

    # --- 2. eigcon：算本征函数并入库 ---
    _feed(
        exe_dir / "2caltoeigen_levels.exe",
        [str(cfg.jcom), str(model_path), str(cfg.dmax_km),
         "res1", "eig1", dbname],
        cwd=work_dir,
    )

    # eigcon 可能留下名为 '-p' 的目录，安全清理
    stray = work_dir / "-p"
    if stray.is_dir():
        shutil.rmtree(stray)

    # --- 3. eigen2asc：导出 ASCII 本征函数 ---
    subprocess.run(
        [str(exe_dir / "3eigentoasc_levels.exe"),
         str(cfg.nmin), str(cfg.nmax), str(cfg.lmin), str(cfg.lmax),
         dbname, dbname],
        cwd=work_dir, check=True, capture_output=True, text=True,
    )
    stray = work_dir / "-p"
    if stray.is_dir():
        shutil.rmtree(stray)

    # res1 复制为可解析文本
    shutil.copy(work_dir / "res1", work_dir / "res1.txt")
    return work_dir / "res1.txt", work_dir / dbname


# ---------------------------------------------------------------- 3. 解析模式表与本征函数


def parse_mode_table(path):
    """稳健解析 mineos 模式表（res1 / res1.txt）。

    自动定位 'mode ... phs vel' 表头；剥离模式序号的 's'/'t' 后缀；
    跳过 '***' 标记行。返回结构数组：
    n, l, phs_vel, freq_mhz, period_s, grp_vel, q, raylquo
    """
    rows = []
    in_table = False
    header_seen = False
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not header_seen and s.startswith("mode") and "phs" in s:
                header_seen = True
                in_table = True
                continue
            # 有表头时只吃表头之后的行；无表头（文件已被预处理）时全文扫描
            if header_seen and not in_table:
                continue
            if not s or s.startswith("*") or s.startswith("mode"):
                continue
            if "integration precision" in s or "gravity cut off" in s:
                continue
            tok = s.split()
            if len(tok) < 8:
                continue
            # 原始 res1 中 n 与 l 之间有一列模式类型 (s/t)
            if tok[1] in ("s", "t", "S", "T"):
                tok = [tok[0]] + tok[2:]
            if len(tok) < 8:
                continue
            try:
                n = int(tok[0].rstrip("stST"))
                rows.append((n, int(tok[1]), float(tok[2]), float(tok[3]),
                             float(tok[4]), float(tok[5]), float(tok[6]),
                             float(tok[7])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"未能在 {path} 中解析到任何模式，请检查文件。")
    dt = np.dtype([("n", int), ("l", int), ("phs", float), ("f_mhz", float),
                   ("period", float), ("grp", float), ("q", float),
                   ("raylquo", float)])
    return np.array(rows, dtype=dt)


def load_spheroidal_eigenfunction(asc_path, l):
    """读取球型模式 ASC 文件，返回表面在末端的 U, V 数组。

    ASC 列：radius, u, up, v, vp, p, pp（半径从伪半径 6371 km 向中心递减，
    此处反转为 中心->表面）。V 乘 sqrt(l(l+1)) 还原 mineos 的 1/sqrt(l(l+1)) 约定。
    """
    arr = np.loadtxt(asc_path, skiprows=1)
    U = arr[::-1, 1]
    V = arr[::-1, 3] * np.sqrt(l * (l + 1))
    return U, V


# ---------------------------------------------------------------- 4. 响应函数计算


class ResponseCalculator:
    """由模式库 + 加密模型计算 l 阶响应函数（逻辑对齐 TransfCalc_2025new）。"""

    def __init__(self, folder, code_file, mode_table="res1.txt",
                 db_dir="db25new", l_degree=2,
                 static_rad=0.5, static_hor=0.25):
        """
        static_rad/static_hor: B 类（引潮）响应到A类的变换系数（乘 R）。
        默认 0.5 / 0.25 为 l=2 的值；|T_B - coef*R| 即动力学 B 类响应，
        其量级与 A 类（Dyson 体力）响应相近，便于对比。
        """
        self.l = l_degree
        self.k = np.sqrt(l_degree * (l_degree + 1))
        self.static_rad = static_rad
        self.static_hor = static_hor
        self.folder = Path(folder)

        model = np.loadtxt(self.folder / code_file)
        self.r, self.rho, self.vs = model[:, 0], model[:, 1], model[:, 3]
        self.mu = self.rho * self.vs ** 2
        self.R = self.r[-1]

        tab = parse_mode_table(self.folder / mode_table)
        self.modes = tab[tab["l"] == l_degree]
        self.db_dir = self.folder / db_dir

        self._precompute()

    def _eig_path(self, n):
        return self.db_dir / f"S.{n:07d}.{self.l:07d}.ASC"

    def _precompute(self):
        r, rho, mu, k, l = self.r, self.rho, self.mu, self.k, self.l
        d_mu = np.gradient(mu, r)
        props = []
        for m in self.modes:
            path = self._eig_path(m["n"])
            if not path.exists():
                continue
            U, V = load_spheroidal_eigenfunction(path, l)
            W = U + (l + 1) * V / k
            I_n = simpson_irregular(rho * (U**2 + V**2) * r**2, r)
            if abs(I_n) < 1e-30:
                continue
            L_tid = simpson_irregular(rho * W * r**3, r)
            L_dys = (simpson_irregular(W * d_mu * r**2, r)
                     - W[-1] * mu[-1] * r[-1] ** 2)
            props.append(dict(
                omega=2 * np.pi * m["f_mhz"] * 1e-3,
                q=max(m["q"], 1.0),
                h=U[-1] * L_tid / I_n,
                kk=(V[-1] / k) * L_tid / I_n,
                A_rad=U[-1] * L_dys / I_n,
                A_hor=(V[-1] / k) * L_dys / I_n,
            ))
        if not props:
            raise RuntimeError("没有可用模式，请检查 ASC 文件与 l 值。")
        self.p = {k_: np.array([d[k_] for d in props])
                  for k_ in ("omega", "q", "h", "kk", "A_rad", "A_hor")}
        self.n_used = len(props)
        self.n_total = len(self.modes)

    def response(self, freq_hz):
        """矢量化叠加所有模式，返回响应字典。"""
        w = 2 * np.pi * np.asarray(freq_hz, dtype=float)
        p = self.p
        denom = p["omega"]**2 - w[:, None]**2 + 1j * w[:, None] * p["omega"] / p["q"]
        TB_r = -np.sum(p["h"] * w[:, None]**2 / (2 * denom), axis=1)
        TB_h = -np.sum(p["kk"] * w[:, None]**2 / (2 * denom), axis=1)
        TA_r = np.sum(p["A_rad"] / denom, axis=1)
        TA_h = np.sum(p["A_hor"] / denom, axis=1)
        return {
            "freq_hz": np.asarray(freq_hz, dtype=float),
            "T_A_rad": np.abs(TA_r),
            "T_A_hor": np.abs(TA_h),
            "T_B_rad": np.abs(TB_r),
            "T_B_hor": np.abs(TB_h),
            "T_B_rad_dyn": np.abs(TB_r - self.static_rad * self.R),
            "T_B_hor_dyn": np.abs(TB_h - self.static_hor * self.R),
        }


def default_freq_grid():
    """默认频率网格：10^-3.5 ~ 10^1 Hz（1000+100 点）。"""
    return 10.0 ** np.concatenate([
        np.linspace(-3.5, -1.5, 1000, endpoint=False),
        np.linspace(-1.5, 1.0, 100),
    ])


def save_response(results, path, l, save_format="compact",
                  static_rad=0.5, static_hor=0.25):
    """保存响应函数。

    save_format:
      "compact" — 5 列（兼容原 RespF 格式）：freq, A类径向, A类切向,
                  B类径向(减系数), B类切向(减系数)
      "full"    — 7 列：额外包含未减系数的 B 类响应 |T_B_r|, |T_B_h|
    两类响应定义：
      A 类（T_A）：Dyson 体力（direct body forcing）驱动下的表面位移响应；
      B 类（T_B）：引潮力（tidal forcing）驱动下的表面位移响应，
                  其变换系数分别为 static_rad*R（径向）和 static_hor*R（切向），
                  减去变换系数后与 A 类量级相近。
    """
    if save_format == "full":
        header = (
            f"Response functions for l={l}\n"
            "Columns:\n"
            "1. Frequency (Hz)\n"
            "2. A-type Radial Response T_A_r (m)   [Dyson body forcing]\n"
            "3. A-type Horizontal Response T_A_h (m)\n"
            "4. B-type Radial Response |T_B_r| (m) [tidal, static NOT subtracted]\n"
            "5. B-type Horizontal Response |T_B_h| (m)\n"
            f"6. B-type Radial, dynamic: |T_B_r - {static_rad}*R| (m)\n"
            f"7. B-type Horizontal, dynamic: |T_B_h - {static_hor}*R| (m)"
        )
        data = np.column_stack([
            results["freq_hz"], results["T_A_rad"], results["T_A_hor"],
            results["T_B_rad"], results["T_B_hor"],
            results["T_B_rad_dyn"], results["T_B_hor_dyn"],
        ])
    else:
        header = (
            f"Response functions for l={l}\n"
            "A-type: Dyson body forcing; B-type: tidal (static limit subtracted)\n"
            "Columns:\n"
            "1. Frequency (Hz)\n"
            "2. A-type Radial Response T_A_r (m)\n"
            "3. A-type Horizontal Response T_A_h (m)\n"
            f"4. B-type Radial |T_B_r - {static_rad}*R| (m)\n"
            f"5. B-type Horizontal |T_B_h - {static_hor}*R| (m)"
        )
        data = np.column_stack([
            results["freq_hz"], results["T_A_rad"], results["T_A_hor"],
            results["T_B_rad_dyn"], results["T_B_hor_dyn"],
        ])
    np.savetxt(path, data, header=header, fmt="%.6e", delimiter="\t")


def plot_response(results, l, prefix, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir = Path(out_dir)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(results["freq_hz"], results["T_A_rad"], "r-",
              label=r"$T^A_r$ (Dyson Forcing)")
    ax.loglog(results["freq_hz"], results["T_B_rad_dyn"], "b--",
              label=r"$|T^B_r - 0.5 R|$")
    ax.set(title=f"Radial Response Functions (l={l})",
           xlabel="Frequency (Hz)", ylabel="Response Amplitude (m)")
    ax.legend(); ax.grid(True, which="both", ls="--")
    fig.savefig(out_dir / f"{prefix}_radial.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(results["freq_hz"], results["T_A_hor"], "r-",
              label=r"$T^A_h$ (Dyson Forcing)")
    ax.loglog(results["freq_hz"], results["T_B_hor_dyn"], "b--",
              label=r"$|T^B_h - 0.25 R|$")
    ax.set(title=f"Horizontal Response Functions (l={l})",
           xlabel="Frequency (Hz)", ylabel="Response Amplitude (m)")
    ax.legend(); ax.grid(True, which="both", ls="--")
    fig.savefig(out_dir / f"{prefix}_horizontal.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- 命令行


def _model_config_from_args(args):
    if args.dr_km is not None:
        return ModelConfig(scheme="dr", dr_km=args.dr_km)
    if args.layers is not None:
        return ModelConfig(scheme="layers", target_layers=args.layers)
    return ModelConfig(scheme="legacy", mult=args.mult, mj=args.mj,
                       crit_dis_km=args.crit)


def cmd_model(args):
    simple = load_simple_model(args.model)
    mcfg = _model_config_from_args(args)
    arr9 = refine_model(simple, mcfg)
    mp, cp = write_mineos_models(arr9, args.name, args.out,
                                 nic=args.nic, noc=args.noc)
    print(f"模型加密完成: {len(arr9)} 层 (方案 {mcfg.scheme}) -> {mp} , {cp}")
    return mp


def cmd_eigen(args):
    ecfg = EigenConfig(jcom=args.jcom, eps=args.eps, wgrav=args.wgrav,
                       lmin=args.lmin, lmax=args.lmax, wmin=args.wmin,
                       wmax=args.wmax, nmin=args.nmin, nmax=args.nmax)
    res1, db = run_eigen_pipeline(args.model_file, ecfg, args.exe_dir,
                                  args.out, dbname=args.db)
    tab = parse_mode_table(res1)
    print(f"本征模式计算完成: {len(tab)} 个模式 -> {res1} , ASC 目录 {db}")


def cmd_response(args):
    calc = ResponseCalculator(args.folder, args.code, mode_table=args.table,
                              db_dir=args.db, l_degree=args.l,
                              static_rad=args.static_rad,
                              static_hor=args.static_hor)
    print(f"有效模式 {calc.n_used}/{calc.n_total}")
    freq = default_freq_grid() if args.freq is None else np.loadtxt(args.freq)
    res = calc.response(freq)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    out_txt = Path(args.out) / args.save
    save_response(res, out_txt, args.l, save_format=args.save_format,
                  static_rad=args.static_rad, static_hor=args.static_hor)
    if not args.no_plot:
        plot_response(res, args.l, Path(args.save).stem, args.out)
    print(f"响应函数已保存: {out_txt} (格式 {args.save_format})")


def cmd_all(args):
    mp = cmd_model(args)
    ecfg = EigenConfig(jcom=args.jcom, eps=args.eps, wgrav=args.wgrav,
                       lmin=args.lmin, lmax=args.lmax, wmin=args.wmin,
                       wmax=args.wmax, nmin=args.nmin, nmax=args.nmax)
    res1, db = run_eigen_pipeline(mp, ecfg, args.exe_dir, args.out,
                                  dbname=args.db)
    tab = parse_mode_table(res1)
    print(f"本征模式: {len(tab)} 个")
    # 将加密模型复制进输出目录，使响应阶段自洽
    code_src = Path(args.out) / f"{args.name}_code.txt"
    args.folder, args.code = args.out, code_src.name
    args.table, args.db = "res1.txt", args.db
    cmd_response(args)


def main(argv=None):
    ap = argparse.ArgumentParser(description="月球响应函数端到端计算")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--out", default=".", help="输出目录")
        # --- 模型插值（三种方案，默认 legacy 兼容旧结果）---
        p.add_argument("--dr-km", type=float, default=None,
                       help="插值方案1: 目标最大层厚 km（如 0.22）")
        p.add_argument("--layers", type=int, default=None,
                       help="插值方案2: 目标总层数（如 8001）")
        p.add_argument("--mult", type=int, default=8000, help="legacy 方案: 厚层份数")
        p.add_argument("--mj", type=int, default=200, help="legacy 方案: 薄层份数")
        p.add_argument("--crit", type=float, default=30.0, help="legacy 方案: 临界层厚 km")
        # --- 流体层（默认自动识别 vs=0 且 rho>0 的最深流体块）---
        p.add_argument("--nic", type=int, default=None,
                       help="流体块下方固体层数（默认自动识别）")
        p.add_argument("--noc", type=int, default=None,
                       help="最外流体层 1-based 层号（默认自动识别）")
        # --- 本征模式 ---
        p.add_argument("--exe-dir", default=None,
                       help="MINEOS exe 目录（默认自动发现 mineos_build/bin）")
        p.add_argument("--db", default="db25new", help="本征函数库目录名")
        p.add_argument("--jcom", type=int, default=3, help="1径向 2环型 3球型 4内核环型")
        p.add_argument("--eps", type=float, default=1e-8, help="积分精度")
        p.add_argument("--wgrav", type=float, default=0.5,
                       help="自引力阈值 mHz（高于此频率忽略重力项）")
        p.add_argument("--lmin", type=int, default=2)
        p.add_argument("--lmax", type=int, default=2)
        p.add_argument("--wmin", type=float, default=0.1, help="频率下限 mHz")
        p.add_argument("--wmax", type=float, default=300.0, help="频率上限 mHz")
        p.add_argument("--nmin", type=int, default=0,
                       help="起始径向阶数（可不从最小 n 开始）")
        p.add_argument("--nmax", type=int, default=400, help="最大径向阶数")
        # --- 响应输出 ---
        p.add_argument("--l", type=int, default=2, help="响应计算的球谐度")
        p.add_argument("--static-rad", type=float, default=0.5,
                       help="B类径向变换系数（乘R，l=2 为 0.5）")
        p.add_argument("--static-hor", type=float, default=0.25,
                       help="B类切向变换系数（乘R，l=2 为 0.25）")
        p.add_argument("--save-format", choices=["compact", "full"],
                       default="compact",
                       help="compact=5列(兼容旧格式) / full=7列(含未减变换系数的B类)")
        p.add_argument("--save", default="response_functions.txt")
        p.add_argument("--freq", default=None, help="自定义频率数组文件")
        p.add_argument("--no-plot", action="store_true")

    p = sub.add_parser("model", help="只做模型加密")
    p.add_argument("--model", required=True, help="简化模型文件")
    p.add_argument("--name", default="RM_out")
    common(p); p.set_defaults(func=cmd_model)

    p = sub.add_parser("eigen", help="只做本征模式计算")
    p.add_argument("--model-file", required=True, help="mineos 卡片模型 .txt")
    common(p); p.set_defaults(func=cmd_eigen)

    p = sub.add_parser("response", help="只做响应函数计算")
    p.add_argument("--folder", required=True, help="模式库文件夹")
    p.add_argument("--code", required=True, help="加密模型数据文件名")
    p.add_argument("--table", default="res1.txt")
    common(p); p.set_defaults(func=cmd_response)

    p = sub.add_parser("all", help="全流程：模型->本征模式->响应曲线")
    p.add_argument("--model", required=True, help="简化模型文件")
    p.add_argument("--name", default="RM_out")
    common(p); p.set_defaults(func=cmd_all)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
