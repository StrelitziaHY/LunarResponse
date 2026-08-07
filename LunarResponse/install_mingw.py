# -*- coding: utf-8 -*-
"""
install_mingw.py — 安装便携版 MinGW-w64 gfortran 工具链（无需管理员权限）

从清华大学 MSYS2 镜像下载 ucrt64 工具链包并解压到 ./tools/msys64。
仅在需要重新编译 MINEOS（见 build_mineos.py）时使用；
如果只运行 lunar_response.py，bin/ 中预编译的 exe 已足够，无需本脚本。

Installs a portable MinGW-w64 gfortran toolchain into ./tools/msys64
(from the Tsinghua MSYS2 mirror; no admin rights required).
Only needed to rebuild MINEOS; the prebuilt binaries in bin/ need no compiler.
"""
import io
import re
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent.resolve()
REPO = "https://mirrors.tuna.tsinghua.edu.cn/msys2/mingw/ucrt64/"
DEST = HERE / "tools" / "msys64"
PKGS = [
    "mingw-w64-ucrt-x86_64-gcc-fortran",
    "mingw-w64-ucrt-x86_64-gcc",
    "mingw-w64-ucrt-x86_64-gcc-libs",
    "mingw-w64-ucrt-x86_64-gcc-libgfortran",
    "mingw-w64-ucrt-x86_64-binutils",
    "mingw-w64-ucrt-x86_64-crt-git",
    "mingw-w64-ucrt-x86_64-headers-git",
    "mingw-w64-ucrt-x86_64-libwinpthread-git",
    "mingw-w64-ucrt-x86_64-winpthreads-git",
    "mingw-w64-ucrt-x86_64-windows-default-manifest",
    "mingw-w64-ucrt-x86_64-gmp",
    "mingw-w64-ucrt-x86_64-mpfr",
    "mingw-w64-ucrt-x86_64-mpc",
    "mingw-w64-ucrt-x86_64-isl",
    "mingw-w64-ucrt-x86_64-zlib",
    "mingw-w64-ucrt-x86_64-zstd",
    "mingw-w64-ucrt-x86_64-libiconv",
]


def ensure_zstandard():
    try:
        import zstandard
        return zstandard
    except ImportError:
        print("安装 zstandard 解压库（pip，清华镜像）...")
        libs = HERE / "tools" / "pylibs"
        subprocess.run([sys.executable, "-m", "pip", "install", "--target",
                        str(libs), "-i",
                        "https://pypi.tuna.tsinghua.edu.cn/simple",
                        "zstandard", "-q"], check=True)
        sys.path.insert(0, str(libs))
        import zstandard
        return zstandard


def resolve_names():
    html = urllib.request.urlopen(REPO, timeout=30).read().decode("utf-8", "replace")
    found = {}
    for stem in PKGS:
        pat = re.compile(re.escape(stem) + r"-[0-9][^\"']*?-any\.pkg\.tar\.zst")
        matches = [m for m in pat.findall(html) if not m.endswith(".sig")]
        if not matches:
            print("!! not found:", stem)
            continue
        found[stem] = sorted(matches)[-1]
    return found


def fetch(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)


def extract_zst(zstd, archive, dest):
    dctx = zstd.ZstdDecompressor()
    with open(archive, "rb") as fh:
        with dctx.stream_reader(fh) as reader:
            data = io.BytesIO(reader.read())
    with tarfile.open(fileobj=data) as tf:
        tf.extractall(dest, filter="data")


def main():
    zstd = ensure_zstandard()
    DEST.mkdir(parents=True, exist_ok=True)
    pkgs = resolve_names()
    print(f"resolved {len(pkgs)}/{len(PKGS)} packages")
    for stem, fname in pkgs.items():
        arc = DEST / fname
        if not arc.exists():
            print("download", fname, flush=True)
            fetch(REPO + fname, arc)
        print("extract ", fname, flush=True)
        extract_zst(zstd, arc, DEST)
        arc.unlink()
    gfortran = DEST / "ucrt64" / "bin" / "gfortran.exe"
    r = subprocess.run([str(gfortran), "--version"], capture_output=True, text=True)
    print(r.stdout.splitlines()[0] if r.returncode == 0 else "!! gfortran 自检失败")
    print("done ->", DEST)


if __name__ == "__main__":
    main()
