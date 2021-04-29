#!/usr/bin/env python
from __future__ import print_function

import argparse
from collections import defaultdict
from configparser import ConfigParser
import os
import re
import sys
import platform


# CPU_ARCHS and GPU_ARCHS are alphabetically ordered
CPU_ARCHS = frozenset(
    (
        "BDW",
        "EPYC",
        "HSW",
        "SKX",
        "WSM",
    ),
)

GPU_ARCHS = frozenset(
    (
        "Kepler30",
        "Maxwell50",
        "Pascal60",
        "Pascal61",
        "Turing75",
        "Volta70",
        "Volta72",
    ),
)

KNOWN_ARCHS = {"CPU": CPU_ARCHS, "GPU": GPU_ARCHS}
DEFAULT_ARCHS = {"CPU": "BDW", "GPU": "Pascal60"}

_GPU_FLAG_DEPRECATION_MESSAGE = (
    "The -gpu flag is deprecated. Using it will raise an error in a future release. "
    "Please explicitly request a GPU architecture via the -arch argument."
)


def get_user_config_file():
    """
    Look for a file named "idefix.cfg" and return the first one found from
    1) the directory of invocation (cwd)
    2) $XDG_CONFIG_HOME (defaults to $HOME/.config)
    return `None` otherwise.
    """
    if os.path.isfile("idefix.cfg"):
        return "idefix.cfg"
    else:
        if platform.system().lower().startswith("win"):
            # Windows
            env_var = "APPDATA"
            default_usr_dir = "AppData"
        else:
            # POSIX
            env_var = "XDG_CONFIG_HOME"
            default_usr_dir = ".config"

        config_dir = os.environ.get(
            env_var,
            os.path.join(os.path.expanduser("~"), default_usr_dir),
        )
        conf_file = os.path.join(config_dir, "idefix.cfg")
        if os.path.isfile(conf_file):
            return conf_file
    return None


def _add_parser_args(parser):
    parser.add_argument(
        "directory", nargs="?", default=os.getcwd(), help="target directory",
    )

    parser.add_argument("-mhd", action="store_true", help="enable MHD")

    parser.add_argument(
        "-gpu",
        action="store_true",
        help=_GPU_FLAG_DEPRECATION_MESSAGE,
    )

    parser.add_argument("-cxx", help="override default compiler")

    parser.add_argument(
        "-arch",
        nargs="+",
        dest="archs",
        choices=CPU_ARCHS.union(GPU_ARCHS),
        help="target Kokkos architecture (accepts up to one CPU and up to one GPU archs)",
    )
    parser.add_argument(
        "-openmp",
        help="enable OpenMP parallelism (not available with GPU architectures)",
        action="store_true",
    )
    parser.add_argument("-mpi", action="store_true", help="enable MPI parallelism")
    parser.add_argument("-defs", help="specify a custom name for definitions.hpp")


def is_gpu_requested(requested_archs):
    return any((a in GPU_ARCHS for a in requested_archs.values()))


def _set_arch_from_conf(config_file_archs, user_config, kind):
    # valid kinds are "CPU" and "GPU"
    known_names = KNOWN_ARCHS[kind]

    name = user_config.get(kind)
    if name is None:
        return

    if name not in known_names:
        print(
            "Warning: unknown {} arch '{}' found in configuration file will be ignored.".format(
                kind,
                name,
            ),
            file=sys.stderr,
        )
        return
    config_file_archs[kind] = name


def parse_archs(cli_archs, user_config):
    """
    Parse user-requested architectures for CPU and GPU, from the command line
    and the persistent user configuration file.
    At most one of each arch kind (CPU, GPU) may be requested.
    CLI arguments (`cli_archs`) take priority of the user config file (`user_config`).

    Returns
    -------
    selected_arch: dict
      keys are "CPU" and "GPU". Both are optional (may be missing).
    """

    if cli_archs is None:
        cli_archs = []
    if len(cli_archs) > 2:
        raise ValueError(
            "Error: received more than two architectures ({}).".format(
                ", ".join(cli_archs),
            ),
        )

    requested_archs = {}

    # parse user configuration (if any)
    if user_config:
        config_file_archs = {}
        _set_arch_from_conf(config_file_archs, user_config, "CPU")
        _set_arch_from_conf(config_file_archs, user_config, "GPU")
        requested_archs.update(config_file_archs)

    # parse command line arguments (take priority over everything else)
    for arch_type, archs in KNOWN_ARCHS.items():
        vals = sorted(list(archs.intersection(set(cli_archs))))
        if not vals:
            continue
        if len(vals) > 1:
            raise ValueError(
                "Error: received more than one {} architecture ({}).".format(
                    arch_type, ", ".join(vals),
                ),
            )

        requested_archs[arch_type] = vals[0]

    # finally, make sure that we always have a CPU architecture since
    # it is used even in GPU setups
    requested_archs.setdefault("CPU", DEFAULT_ARCHS["CPU"])
    return requested_archs


def _get_sed():
    # Build a sed command which is compatible with the current platofm (BSD and GNU diverge on that)
    sed = ""
    if platform.system() == 'Darwin':
        sed = "sed -i '' "
    else:
        sed = "sed -i"

    return sed


def _get_makefile_options(
    archs,
    cxx,
    openmp,
    mpi,
    mhd,
    sed,
    defs,
):
    # using a default dict to allow setting key value pairs as
    # >>> options[key] += value

    options = defaultdict(str)
    options["cxxflags"] = "-O3"
    options["sed-command"] = sed

    if is_gpu_requested(archs):
        options["extraLine"] = '\nKOKKOS_CUDA_OPTIONS = "enable_lambda"'
        options["cxx"] = "${KOKKOS_PATH}/bin/nvcc_wrapper"
        options["kokkosDevices"] = '"Cuda"'
        options["kokkosArch"] = "{CPU},{GPU}".format(**archs)

        # Enforce backend compiler for nvcc
        nvcc = "\nexport NVCC_WRAPPER_DEFAULT_COMPILER = {}"
        if cxx:
            options["extraLine"] += nvcc.format(cxx)
        elif mpi:
            options["extraLine"] += nvcc.format("mpicxx")
    else:
        if cxx:
            options["cxx"] = cxx
        elif mpi:
            options["cxx"] = "mpicxx"
        else:
            options["cxx"] = "g++"

        options["kokkosArch"] = archs["CPU"]
        options["kokkosDevices"] = '"OpenMP"' if openmp else '"Serial"'

    if mpi:
        options["extraIncludeDir"] += " -I$(SRC)/dataBlock/mpi"
        options["extraVpath"] += ":$(SRC)/dataBlock/mpi"
        options["extraObj"] += " mpi.o"
        options["cxxflags"] += " -DWITH_MPI"

    if mhd:
        options["extraIncludeDir"] += " -I$(SRC)/hydro/MHDsolvers"
        options["extraVpath"] += ":$(SRC)/hydro/MHDsolvers"
        options["cxxflags"] += " -DMHD=YES"
    else:
        options["extraIncludeDir"] += " -I$(SRC)/hydro/HDsolvers"
        options["extraVpath"] += ":$(SRC)/hydro/HDsolvers"
        options["cxxflags"] += " -DMHD=NO"

    if defs:
        options["cxxflags"] += " -DDEFINITIONS_FILE=\'\"" + defs + "\"\'"

    return options


def _write_makefile(
    directory,
    options,
):
    with open(os.path.join(os.environ["IDEFIX_DIR"], "Makefile.in")) as fh:
        data = fh.read()

    # apply subsitutions
    for key, val in options.items():
        data = data.replace(r"@{}@".format(key), val)

    # cleanup unused place holders
    data = re.sub(r"@.+@", "", data)

    with open(os.path.join(directory, "Makefile"), "w") as fh:
        fh.write(data)


def _get_report(
    archs,
    openmp,
    mpi,
    mhd,
    makefile_options,
):
    def status(name, flag):
        prefix = "en" if flag else "dis"
        return "{}: {}abled".format(name, prefix)

    report_lines = []
    report_lines += [
        "-----------------------------------------------------------",
        "Idefix succesfully configured with the following options:",
        "",
        status("MHD", mhd),
        "Compiler: {}".format(makefile_options["cxx"]),
        status("MPI", mpi),
    ]

    arch_type = "GPU" if is_gpu_requested(archs) else "CPU"
    report_lines += [
        "Execution target: {}".format(arch_type),
        "Target architecture: {}".format(archs[arch_type]),
    ]
    if arch_type == "CPU":
        report_lines.append(status("OpenMP", openmp))

    report_lines += [
        "Cflags: {}".format(makefile_options["cxxflags"]),
        "-----------------------------------------------------------",
    ]

    return "\n".join(report_lines)


def main(argv=None):

    if "IDEFIX_DIR" not in os.environ:
        print(
            "Error: IDEFIX_DIR environment variable is not defined.",
            file=sys.stderr,
        )
        return 1

    # get the persistent configuration options
    user_config = {}
    config_file = get_user_config_file()
    if config_file is not None:
        cf = ConfigParser()
        try:
            cf.read(config_file)
            # we only care about one relevant section from the config file
            user_config = cf["compilation"]
        except KeyError:
            print(
                (
                    "Warning: located a configuration file '{}', ".format(config_file),
                    "but couldn't find a 'compilation' section therein.",
                ),
                file=sys.stderr,
            )

    parser = argparse.ArgumentParser("setup")
    _add_parser_args(parser)

    args = parser.parse_args(argv)

    try:
        requested_archs = parse_archs(cli_archs=args.archs, user_config=user_config)
    except ValueError as err:
        print(err, file=sys.stderr)
        return 1

    if args.gpu:
        print("Warning: " + _GPU_FLAG_DEPRECATION_MESSAGE, file=sys.stderr)
        if not is_gpu_requested(requested_archs):
            print(
                "Warning: -gpu flag was received, but no GPU architecture was specified. "
                "Defaulting to %s" % DEFAULT_ARCHS["GPU"],
                file=sys.stderr,
            )
            requested_archs.setdefault("GPU", DEFAULT_ARCHS["GPU"])

    if args.openmp and is_gpu_requested(requested_archs):
        print("Warning: with a GPU arch, -openmp flag is ignored.", file=sys.stderr)
        args.openmp = False

    mysed = _get_sed()

    makefile_options = _get_makefile_options(
        archs=requested_archs,
        cxx=args.cxx or user_config.get("CXX"),
        openmp=args.openmp,
        mpi=args.mpi,
        mhd=args.mhd,
        sed=mysed,
        defs=args.defs,
    )
    try:
        _write_makefile(args.directory, makefile_options)
    except OSError:
        filename = os.path.join(args.directory, "Makefile")
        print("Error: could not write to {}".format(filename), file=sys.stderr)
        return 1

    report = _get_report(
        archs=requested_archs,
        openmp=args.openmp,
        mpi=args.mpi,
        mhd=args.mhd,
        makefile_options=makefile_options,
    )
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
