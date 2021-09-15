# this can be run from the top level of the repo with
# python3 -m pytest test/test_configure.py
# see https://docs.pytest.org/en/stable/usage.html#cmdline
import os
from pathlib import Path
import textwrap
import pytest

from configure import main, GPU_ARCHS, CPU_ARCHS  # noqa: E402

IDEFIX_DIR = str(Path(__file__).parents[1])


@pytest.fixture()
def in_tmp_dir(tmp_path):
    previous_dir = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield
    finally:
        os.chdir(previous_dir)


@pytest.fixture()
def no_config_file(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_main_no_idefix(capsys, monkeypatch, tmp_path):
    monkeypatch.delenv("IDEFIX_DIR", raising=False)
    ret = main([])
    assert ret != 0

    out, err = capsys.readouterr()
    assert out == ""
    assert err == "Error: IDEFIX_DIR environment variable is not defined.\n"
    assert not Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_parse_archs_overload(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", "BDW", "Kepler30", "HSW"])
    assert ret != 0

    out, err = capsys.readouterr()
    assert out == ""
    assert err == "Error: received more than two architectures (BDW, Kepler30, HSW).\n"
    assert not Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_parse_archs_two_cpus(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", "BDW", "HSW"])
    assert ret != 0

    out, err = capsys.readouterr()
    assert out == ""
    assert err == "Error: received more than one CPU architecture (BDW, HSW).\n"
    assert not Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_parse_archs_two_gpus(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", "Kepler30", "Maxwell50"])
    assert ret != 0

    out, err = capsys.readouterr()
    assert out == ""
    assert err == "Error: received more than one GPU architecture (Kepler30, Maxwell50).\n"
    assert not Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_main_default_success(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main([])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Idefix succesfully configured" in out
    assert "Execution target: CPU" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.parametrize("arch", GPU_ARCHS)
@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_main_auto_gpu_mode_solo(arch, capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", arch])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.parametrize("arch", CPU_ARCHS)
@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_main_auto_cpu_mode_solo(arch, capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", arch])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: CPU" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "no_config_file")
def test_main_gpu_flag_deprecation(capsys, monkeypatch):
    from configure import _GPU_FLAG_DEPRECATION_MESSAGE
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)
    ret = main(["-arch", "Volta70", "-gpu"])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU" in out
    assert err == "Warning: %s\n" % _GPU_FLAG_DEPRECATION_MESSAGE
    assert Path("Makefile").is_file()


@pytest.fixture()
def tmp_config_file_gpu(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    file = tmp_path / "idefix.cfg"
    with open(file, "w") as fh:
        fh.write(
            textwrap.dedent(
                """
                [compilation]
                GPU = Turing75
                """,
            ),
        )
    return file


@pytest.mark.usefixtures("in_tmp_dir", "tmp_config_file_gpu")
def test_main_gpu_from_conf_file(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)

    # explicitly request a CPU arch on the command line
    ret = main(["-arch", "EPYC"])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU\n" in out
    assert "Target architecture: Turing75\n" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "tmp_config_file_gpu")
def test_main_arch_from_conf_file(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)

    ret = main([])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU\n" in out
    assert "Target architecture: Turing75\n" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "tmp_config_file_gpu")
def test_main_arch_cli_override_conf_file(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)

    ret = main(["-arch", "Volta70"])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU\n" in out
    assert "Target architecture: Volta70\n" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir", "tmp_config_file_gpu")
def test_main_arch_local_override_global_conf_file(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)

    with open("idefix.cfg", "w") as fh:
        fh.write(
            textwrap.dedent(
                """
                [compilation]
                GPU = Pascal60
                """,
            ),
        )

    ret = main([])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Execution target: GPU\n" in out
    assert "Target architecture: Pascal60\n" in out
    assert err == ""
    assert Path("Makefile").is_file()


@pytest.mark.usefixtures("in_tmp_dir")
def test_main_set_compiler_from_conf_file(capsys, monkeypatch):
    monkeypatch.setenv("IDEFIX_DIR", IDEFIX_DIR)

    with open("idefix.cfg", "w") as fh:
        fh.write(
            textwrap.dedent(
                """
                [compilation]
                CXX = MY_COMPILER
                """,
            ),
        )

    ret = main([])
    assert ret == 0

    out, err = capsys.readouterr()
    assert "Compiler: MY_COMPILER\n" in out
    assert err == ""
    assert Path("Makefile").is_file()
