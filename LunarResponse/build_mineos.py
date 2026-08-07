# -*- coding: utf-8 -*-
"""
build_mineos.py — MINEOS 一键构建脚本（含层数上限 mk 调整）

用法：
    python build_mineos.py                 # 默认 mk=50000，编译到 ./bin
    python build_mineos.py --mk 100000     # 调整层数上限，一条命令完成
    python build_mineos.py --mk 80000 --out bin_mk80000

原理：
    src/ 中的三个 .f 文件里所有 `parameter (mk=N)`（共 17 处，必须保持一致）
    被统一替换为目标值后写入临时构建目录编译，src/ 原始源码永不修改。
    默认尝试完全静态链接（-static-libgcc -static-libgfortran -static-libquadmath），
    成功后 exe 不依赖任何 MinGW DLL，可单独分发。
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent.resolve()
SRC = HERE / "src"
PROGRAMS = [  # (输出名, [源文件])
    ("1calc_levels.exe", ["minos_bran.f"]),
    ("2caltoeigen_levels.exe", ["eigcon.f"]),
    ("3eigentoasc_levels.exe", ["eigen2asc.f", "fdb_eigen.f"]),
]
DLLS = ["libgfortran-5.dll", "libquadmath-0.dll", "libgcc_s_seh-1.dll"]
MK_PAT = re.compile(r"parameter\s*\(\s*mk\s*=\s*\d+\s*\)", re.IGNORECASE)


def find_gfortran():
    """按 PATH -> 本地工具链顺序找 gfortran。"""
    if shutil.which("gfortran"):
        return "gfortran", None
    for root in (HERE, HERE.parent):  # 包内 tools/ 或上级 tools/
        local = root / "tools" / "msys64" / "ucrt64" / "bin"
        if (local / "gfortran.exe").exists():
            return str(local / "gfortran.exe"), local
    raise SystemExit("未找到 gfortran。请先运行 install_mingw.py 或把 gfortran 加入 PATH。")


def patch_mk(src_text, mk):
    """替换全部 parameter (mk=...) 为指定值。"""
    new, n = MK_PAT.subn(f"parameter (mk={mk})", src_text)
    return new, n


def exe_imports(exe_path):
    """扫描 exe 的 DLL 导入表。"""
    data = Path(exe_path).read_bytes()
    return sorted({m.decode() for m in
                   re.findall(rb"[A-Za-z0-9_\-]+\.[Dd][Ll][Ll]", data)})


def main():
    ap = argparse.ArgumentParser(description="MINEOS 构建（mk 可调）")
    ap.add_argument("--mk", type=int, default=50000, help="层数上限（默认 50000）")
    ap.add_argument("--out", default=str(HERE / "bin"), help="exe 输出目录")
    ap.add_argument("--no-static", action="store_true", help="禁用静态链接尝试")
    args = ap.parse_args()

    gfortran, local_bin = find_gfortran()
    out_dir = Path(args.out).resolve()  # 相对路径须在调用者的 cwd 下解析
    out_dir.mkdir(parents=True, exist_ok=True)

    # gfortran 驱动依赖同目录工具/DLL：本地工具链须注入 PATH
    env = None
    if local_bin is not None:
        import os
        env = os.environ.copy()
        env["PATH"] = str(local_bin) + os.pathsep + env.get("PATH", "")

    with tempfile.TemporaryDirectory(dir=HERE) as tmp:
        tmp = Path(tmp)
        shutil.copytree(SRC / "fdb", tmp / "fdb")
        shutil.copy(SRC / "fdb_eigen.h", tmp / "fdb_eigen.h")  # fdb_eigen.f 不带目录前缀
        n_patched = 0
        for f in ("minos_bran.f", "eigcon.f", "eigen2asc.f", "fdb_eigen.f"):
            text = (SRC / f).read_text(encoding="utf-8", errors="replace")
            text, n = patch_mk(text, args.mk)
            n_patched += n
            (tmp / f).write_text(text, encoding="utf-8")
        print(f"mk={args.mk}: 共替换 {n_patched} 处 parameter (mk=...)")

        flags = ["-O2"]
        if not args.no_static:
            # 注：本 MSYS2 工具链的静态 libgfortran.a 链接即崩溃，
            # 故 libgfortran 保持动态（构建后自动附带 libgfortran-5.dll）。
            flags += ["-static-libgcc", "-static-libquadmath"]
        for exe, srcs in PROGRAMS:
            cmd = [gfortran, *flags, "-o", str(out_dir / exe), *srcs]
            r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, env=env)
            if r.returncode != 0:
                print("诊断信息：")
                print("  gfortran =", gfortran)
                print("  cmd =", cmd)
                print("  returncode =", r.returncode)
                print("  stdout =", repr(r.stdout[-1500:]))
                print("  stderr =", repr(r.stderr[-1500:]))
                raise SystemExit(f"编译失败: {exe}")
            print("编译成功:", exe)

    # 运行时 DLL 无条件全部附带：libgfortran-5.dll 自身还依赖 libquadmath/libgcc，
    # 仅按 exe 导入表复制会漏掉二级依赖导致 0xC0000135。
    dll_dir = local_bin or Path(shutil.which("gfortran")).parent
    for d in DLLS:
        src_dll = dll_dir / d
        if src_dll.exists() and not (out_dir / d).exists():
            shutil.copy(src_dll, out_dir / d)
    for exe, _ in PROGRAMS:
        imports = exe_imports(out_dir / exe)
        direct = [i for i in imports if i.lower().startswith("lib")]
        print(f"  {exe}: 直接依赖 {direct}（三个运行时 DLL 已全部附带）")

    # 冒烟测试：程序应能启动并提示输入模型文件名
    p = subprocess.run([str(out_dir / "1calc_levels.exe")], input="\n",
                       capture_output=True, text=True, timeout=30)
    if "model file" in (p.stdout + p.stderr):
        print("冒烟测试通过：1calc_levels.exe 正常启动。")
    else:
        print("!! 冒烟测试异常，输出：", (p.stdout + p.stderr)[:300])
    print(f"\n完成。exe 位于 {out_dir}（层数上限 mk={args.mk}）")


if __name__ == "__main__":
    main()
