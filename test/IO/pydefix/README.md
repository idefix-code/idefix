This directory shows how to use Idefix with pure python initial conditions and outputs. It reproduces the 2D OrszagTang vortex available in MHD/OrszagTang without requiring any single line of C++ code from the user.

The python script `pydefix_example.py` initializes the initial condition of the OT test (`initflow`) and produces a series of PNG files through matplotlib (`output`).

# Python modules installation

In order to use pydefix, you need to be working in a python environement that includes pybind11. The simplest way is to install the suggested list of packages (that you can deploy in dedicated venv environement)

```bash
pip install -r python_requirements.txt
```

# Configuration

In order to use Pydefix, you need to switch the `Idefix_PYTHON` to ON in cmake. In this particular setup, it is automatically done for you in  CMakeLists.txt to avoid mistakes.

# Running

Just run idefix as usual.

# Troubleshooting

It during configuration stage, you get:

`
CMake Error at CMakeLists.txt:122 (find_package):
  By not providing "Findpybind11.cmake" in CMAKE_MODULE_PATH this project has
  asked CMake to find a package configuration file provided by "pybind11",
  but CMake did not find one.`

It means that cmake cannot find the location of pybind11 (this typically happens on MacOs). In order to locate pybind11, open a python interpreter, and get pybind11 install dir through:

```python
import pybind11
print(pybind11.__file__)
```

You can then exit the interpreter and set the pybind11_DIR environement variable to the right path:

```bash
export pybind11_DIR=env/lib/python3.10/site-packages/pybind11
```

you can then run cmake which should be able to find pybind11, and compile the code.
