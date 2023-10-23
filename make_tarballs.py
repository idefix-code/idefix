import argparse
import os
import subprocess
import sys
import tarfile
from contextlib import contextmanager
from copy import copy
from tempfile import TemporaryDirectory

IDEFIX_DIR = os.path.dirname(__file__)

CORE_EXCLUDE_LIST = [
    ".*",
    "CPPLINT.cfg",
    "make_tarballs.py",
]

KOKKOS_EXCLUDE_LIST = [
    ".*",
    "appveyor.yml",
    "scripts",
]


@contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(previous_dir)


def get_idefix_version_string():
    return (
        subprocess.check_output(["git", "describe", "--tags", "--always"])
        .decode()
        .strip()
    )


def _make_git_archive(source_dir, *, output_tarfile, exclude_list=None, prefix=None):
    cmd = ["git", "archive", "HEAD", "--output", output_tarfile]

    if exclude_list is not None:
        cmd.extend([":(exclude)%s" % _ for _ in exclude_list])

    if prefix is not None:
        cmd.extend(["--prefix", prefix])

    with pushd(source_dir):
        subprocess.check_call(cmd)


def _make_self_contained_tarball(output_dir, *, suffix="", exclude_list=None):

    tarball = "idefix_%s%s.tar.gz" % (get_idefix_version_string(), suffix)
    final_tarfile = os.path.join(output_dir, tarball)
    print("Creating %s ... " % final_tarfile, end="", flush=True)

    idefix_exclude_list = copy(CORE_EXCLUDE_LIST)
    if exclude_list is not None:
        idefix_exclude_list.extend(exclude_list)

    tasks = (
        ("idefix", {"source_dir": IDEFIX_DIR, "exclude_list": idefix_exclude_list}),
        (
            "kokkos",
            {
                "source_dir": os.path.join(IDEFIX_DIR, "src", "kokkos"),
                "exclude_list": KOKKOS_EXCLUDE_LIST,
                "prefix": os.path.join("src", "kokkos") + os.path.sep,
            },
        ),
    )
    with TemporaryDirectory() as work_dir:
        for name, params in tasks:
            output_tarfile = os.path.join(work_dir, "%s.tar" % name)
            _make_git_archive(**params, output_tarfile=output_tarfile)
            with tarfile.open(output_tarfile) as fh:
                fh.extractall(work_dir)
            os.remove(output_tarfile)

        with pushd(work_dir):
            with tarfile.open(tarball, "w:gz") as fh:
                fh.add(".")

        os.replace(os.path.join(work_dir, tarball), final_tarfile)

    filesize = os.path.getsize(final_tarfile)
    print("Done ! Final file size is %.1f MB" % (filesize / 1024 ** 2))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", default=".", nargs="?")
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help=(
            "allow to build tarballs with uncommited changes. "
            "Useful to test this script"
        ),
    )
    args = parser.parse_args(argv)

    with pushd(IDEFIX_DIR):
        # ensure kokkos version is clean
        subprocess.check_call(["git", "submodule", "update", "--init"])

        # break if local copy is dirty
        stdout = subprocess.check_output(["git", "diff", "--name-only"]).decode()
        if not args.allow_dirty and stdout.strip() != "":
            print(
                "ERROR: Local copy is dirty. Commit or stash changes before publishing.",
                file=sys.stderr,
            )
            return 1

    _make_self_contained_tarball(args.output_dir)
    _make_self_contained_tarball(args.output_dir, suffix="_lite", exclude_list=["test"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
