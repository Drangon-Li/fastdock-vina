#!/usr/bin/env python3
"""Build the Vina command-line tools with GCC or Clang and installed Boost."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import shlex
import subprocess


ROOT = Path(__file__).resolve().parents[1]
LIB_SOURCES = "ad4cache cache non_cache conf_independent coords grid szv_grid model monte_carlo mutate parallel_mc parse_pdbqt quasi_newton quaternion random utils vina".split()
LIBRARIES = "system thread serialization filesystem program_options".split()
FLAGS = ["-O3", "-DNDEBUG", "-std=c++14", "-pthread", "-fvisibility=hidden",
         '-DVERSION="v1.2.7-27-g3c65c0b"']


def run(command):
    print(shlex.join(command), flush=True)
    subprocess.run(command, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boost-prefix", type=Path,
                        help="Boost prefix containing include/ and lib/ (e.g. brew --prefix boost); omit for system Boost")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"),
                        help="compiler executable (default: CXX environment variable or c++)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="number of concurrent compiler processes (default: 4)")
    parser.add_argument("--output-dir", type=Path, default=Path("build"))
    parser.add_argument("--with-split", action="store_true",
                        help="also build vina_split")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    output = args.output_dir.resolve()
    flags = FLAGS + ["-I" + str(ROOT / "src/lib")]
    if args.boost_prefix:
        prefix = args.boost_prefix.resolve()
        if not (prefix / "include/boost/version.hpp").is_file():
            parser.error("--boost-prefix must contain include/boost/version.hpp")
        flags += ["-I" + str(prefix / "include")]

    objects_dir = output / "obj"
    objects_dir.mkdir(parents=True, exist_ok=True)
    sources = [ROOT / "src/main/main.cpp"]
    sources += [ROOT / "src/lib" / (name + ".cpp") for name in LIB_SOURCES]
    if args.with_split:
        sources.append(ROOT / "src/split/split.cpp")

    def compile_one(source):
        obj = objects_dir / (source.stem + ".o")
        run([args.cxx] + flags + ["-c", str(source), "-o", str(obj)])
        return str(obj)

    try:
        probe = [args.cxx] + flags + ["-dM", "-E", "-x", "c++", "-"]
        macros = subprocess.check_output(probe, input="#include <boost/version.hpp>\n", text=True)
        boost_version = next((int(line.split()[2]) for line in macros.splitlines()
                              if line.startswith("#define BOOST_VERSION ")), 0)
        if not boost_version:
            raise ValueError("could not read BOOST_VERSION from boost/version.hpp")
        # Boost 1.89 removed the compiled stub for header-only Boost.System.
        names = [name for name in LIBRARIES if name != "system" or boost_version < 108900]
        libraries = ["-lboost_" + name for name in names]
        if args.boost_prefix:
            libdir = prefix / "lib"
            archives = [libdir / ("libboost_" + name + ".a") for name in names]
            if all(path.is_file() for path in archives):
                libraries = [str(path) for path in archives]
            else:
                libraries = ["-L" + str(libdir), "-Wl,-rpath," + str(libdir)] + libraries
        # Rebuild all objects so changed headers, flags and dependencies are honored.
        with ThreadPoolExecutor(max_workers=min(args.jobs, len(sources))) as executor:
            objects = list(executor.map(compile_one, sources))
        vina_objects = objects[:-1] if args.with_split else objects
        run([args.cxx] + flags + vina_objects + libraries + ["-o", str(output / "vina")])
        if args.with_split:
            run([args.cxx] + flags + [objects[-1]] + libraries + ["-o", str(output / "vina_split")])
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, "build failed: " + str(error) + "\n")
    print("Built " + str(output / "vina"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
