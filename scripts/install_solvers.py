"""Download and build/install Concorde and LKH for CI.

Neither solver ships prebuilt for every platform this project is tested
on, and Concorde has no supported native Windows build at all (its own
documentation lists Windows as unsupported without a Cygwin environment),
so the two solvers are handled differently per OS:

- LKH: built from source on Linux/macOS ('make'); the official prebuilt
  'LKH-3.exe' is downloaded directly on Windows.
- Concorde: built from source on Linux/macOS, linked against a prebuilt
  QSopt LP library fetched for the current OS/architecture. Skipped on
  Windows.

Installed executables are copied into '<directory>/bin' and that
directory is prepended to PATH (via $GITHUB_PATH when running in GitHub
Actions, otherwise the caller must add it to PATH itself).
"""

import argparse
import gzip
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
from pathlib import Path


LKH_SOURCE_URL = "http://webhotel4.ruc.dk/~keld/research/LKH-3/LKH-3.0.14.tgz"
LKH_WINDOWS_EXE_URL = "http://webhotel4.ruc.dk/~keld/research/LKH-3/LKH-3.exe"

CONCORDE_SOURCE_URL = "https://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz"
# QSopt static library + header, matched to (platform.system(), platform.machine()).
QSOPT_URLS = {
        ("Linux", "x86_64"): (
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.a",
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.h"),
        ("Darwin", "arm64"): (
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/m1/qsopt.a",
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/m1/qsopt.h"),
        ("Darwin", "x86_64"): (
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/mac64/qsopt.a.gz",
                "https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/mac64/qsopt.h"),
        }


def download(url, path):
    print("Downloading \"" + url + "\" -> \"" + str(path) + "\"")
    import urllib.request
    urllib.request.urlretrieve(url, path)


def run(command, cwd=None):
    print("$ " + " ".join(str(c) for c in command) + (" (cwd=" + str(cwd) + ")" if cwd else ""))
    subprocess.run(command, cwd=cwd, check=True)


def make_executable(path):
    st = os.stat(path)
    os.chmod(path, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def add_to_path(directory):
    github_path = os.environ.get("GITHUB_PATH")
    if github_path:
        with open(github_path, "a") as f:
            f.write(str(directory) + os.linesep)
    else:
        print("(not running in GitHub Actions: add \"" + str(directory) + "\" to PATH yourself)")


def install_lkh(work_dir, bin_dir):
    if platform.system() == "Windows":
        exe_path = bin_dir / "LKH.exe"
        download(LKH_WINDOWS_EXE_URL, exe_path)
    else:
        archive_path = work_dir / "LKH-3.0.14.tgz"
        download(LKH_SOURCE_URL, archive_path)
        with tarfile.open(archive_path) as tar:
            tar.extractall(work_dir)
        source_dir = work_dir / "LKH-3.0.14"
        run(["make"], cwd=source_dir)
        shutil.copy(source_dir / "LKH", bin_dir / "LKH")
        make_executable(bin_dir / "LKH")

    print("LKH installed in \"" + str(bin_dir) + "\".")


def install_qsopt(qsopt_dir, qsopt_a_url, qsopt_h_url):
    qsopt_dir.mkdir(parents=True, exist_ok=True)
    downloaded_a_path = qsopt_dir / Path(qsopt_a_url).name
    download(qsopt_a_url, downloaded_a_path)
    if downloaded_a_path.suffix == ".gz":
        with gzip.open(downloaded_a_path, "rb") as f_in, open(qsopt_dir / "qsopt.a", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        downloaded_a_path.unlink()
    elif downloaded_a_path.name != "qsopt.a":
        downloaded_a_path.rename(qsopt_dir / "qsopt.a")
    download(qsopt_h_url, qsopt_dir / "qsopt.h")


def install_concorde(work_dir, bin_dir):
    system = platform.system()
    if system == "Windows":
        print(
                "Concorde has no supported native Windows build (its "
                "documentation lists Windows as unsupported without a "
                "Cygwin environment); skipping Concorde on Windows.")
        return False

    key = (system, platform.machine())
    if key not in QSOPT_URLS:
        print("No known QSopt build for platform " + str(key) + "; skipping Concorde.")
        return False
    qsopt_a_url, qsopt_h_url = QSOPT_URLS[key]

    qsopt_dir = work_dir / "qsopt"
    install_qsopt(qsopt_dir, qsopt_a_url, qsopt_h_url)

    archive_path = work_dir / "co031219.tgz"
    download(CONCORDE_SOURCE_URL, archive_path)
    with tarfile.open(archive_path) as tar:
        tar.extractall(work_dir)
    source_dir = work_dir / "concorde"
    make_executable(source_dir / "configure")

    configure_command = [
            str(source_dir / "configure"),
            "--with-qsopt=" + str(qsopt_dir.resolve()),
            ]
    if system == "Darwin":
        configure_command.append("--host=darwin")
    run(configure_command, cwd=source_dir)
    run(["make"], cwd=source_dir)

    shutil.copy(source_dir / "TSP" / "concorde", bin_dir / "concorde")
    make_executable(bin_dir / "concorde")

    print("Concorde installed in \"" + str(bin_dir) + "\".")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
            "--directory",
            type=str,
            default="tools",
            help="directory to download and build the solvers in")
    parser.add_argument(
            "--solvers",
            type=str,
            nargs="*",
            default=["concorde", "lkh"],
            help="which solvers to install")
    args = parser.parse_args()

    work_dir = Path(args.directory).resolve()
    bin_dir = work_dir / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)

    installed_any = False
    if "lkh" in args.solvers:
        install_lkh(work_dir, bin_dir)
        installed_any = True
    if "concorde" in args.solvers:
        if install_concorde(work_dir, bin_dir):
            installed_any = True

    if installed_any:
        add_to_path(bin_dir)


if __name__ == "__main__":
    main()
