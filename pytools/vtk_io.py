"""
Created on Mon Nov  3 15:23:00 2014

@author: glesur
"""

import warnings
import numpy as np
import os
import re

# restrict what's included with `import *` to public API
__all__ = [
    "readVTK",
    "readVTKCart",
    "readVTKPolar",
    "readVTKSpherical",
]

# datatype we read
dt = np.dtype(">f")  # Big endian single precision floats
dint = np.dtype(">i4")  # Big endian integer

NATIVE_COORDINATE_REGEXP = re.compile(r"X(1|2|3)(L|C)_NATIVE_COORDINATES")

KNOWN_GEOMETRIES = {
    0: "cartesian",
    1: "polar",
    2: "spherical",
    3: "cylindrical",
}


class VTKDataset(object):
    def __init__(self, filename, geometry=None):
        self.filename = os.path.abspath(filename)
        self.data = {}
        self.native_coordinates = {}
        with open(filename, "rb") as fh:
            self._load_header(fh, geometry=geometry)
            self._load(fh)

        if self.native_coordinates:
            self._setup_coordinates_from_native()

    def _load(self, fh):
        if self._dataset_type in ("RECTILINEAR_GRID", "STRUCTURED_GRID"):
            self._load_hydro(fh)
        elif self._dataset_type == "POLYDATA":
            self._load_particles(fh)
        else:
            raise ValueError("Unknown dataset type '%s'" % self._dataset_type)

    def _load_header(self, fh, geometry=None):
        fh.seek(0)
        # skip over the first 3 lines which normally contains
        # VTK DataFile Version x.x
        # <Comments>
        # BINARY
        for _ in range(3):
            fh.readline()

        line = fh.readline().decode("utf-8")
        assert line.startswith("DATASET ")
        self._dataset_type = line.split()[1]

        self.geometry = geometry

        ref_position = fh.tell()
        line = fh.readline().decode("utf-8")  # DIMENSIONS NX NY NZ or FIELD
        if line.startswith("DIMENSIONS"):
            # Idefix < 0.8
            fh.seek(ref_position)
        elif line.startswith("FIELD"):
            # Idefix >= 0.8
            nfield = int(line.split()[2])
            for _ in range(nfield):
                d = fh.readline().decode("utf-8")
                if d.startswith("GEOMETRY"):
                    geom_flag = np.fromfile(fh, dint, 1)[0]
                    geometry_from_data = KNOWN_GEOMETRIES.get(geom_flag)
                    if geometry_from_data is None:
                        raise RuntimeError(
                            "Found unknown geometry value %d" % geom_flag
                        )
                    elif geometry is not None and geometry_from_data != geometry:
                        warnings.warn(
                            "Received geometry argument '%s' inconsistent with data geometry '%s'."
                            % (geometry, geometry_from_data),
                            stacklevel=3,
                        )
                    self.geometry = geometry_from_data
                elif d.startswith("TIME"):
                    self.t = np.fromfile(fh, dt, 1)
                elif d.startswith("PERIODICITY"):
                    self.periodicity = np.fromfile(fh, dtype=dint, count=3).astype(bool)
                elif NATIVE_COORDINATE_REGEXP.match(d):
                    native_name, _ncomp, native_dim, _dtype = d.split()
                    self.native_coordinates[native_name] = np.fromfile(fh, dtype=dt, count=int(native_dim))
                else:
                    warnings.warn("Found unknown field %s" % d)
                fh.readline()  # skip extra linefeed (empty line)

        if self.geometry is None:
            raise ValueError(
                "No GEOMETRY metadata was found, "
                "please provide a geometry argument to this function"
            )

    def _load_hydro(self, fh):
        s = fh.readline()
        slist = s.decode("utf-8").split()  # DIMENSIONS....
        self.nx = int(slist[1])
        self.ny = int(slist[2])
        self.nz = int(slist[3])

        if self.geometry in ("cartesian", "cylindrical"):
            # CARTESIAN geometry
            # NOTE: cylindrical geometry is meant to be only used in 2D
            #       so the expected coordinates (R, z) are never curvilinear,
            #       which means we can treat them as cartesian
            s = fh.readline()  # X_COORDINATES NX float
            x = np.fromfile(fh, dt, self.nx)

            s = fh.readline()  # Extra line feed added by idefix

            s = fh.readline()  # X_COORDINATES NX float

            y = np.fromfile(fh, dt, self.ny)
            s = fh.readline()  # Extra line feed added by idefix

            s = fh.readline()  # X_COORDINATES NX float

            z = np.fromfile(fh, dt, self.nz)
            s = fh.readline()  # Extra line feed added by idefix
            s = fh.readline()  # POINT_DATA NXNYNZ

            slist = s.split()
            point_type = str(slist[0].decode("utf-8"))
            npoints = int(slist[1])
            s = fh.readline()  # EXTRA LINE FEED

            if point_type == "CELL_DATA":
                # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
                if self.nx > 1:
                    self.nx = self.nx - 1
                    self.x = 0.5 * (x[1:] + x[:-1])
                    # left side of the cell
                    self.xl = x
                else:
                    self.x = x
                    self.xl = x
                if self.ny > 1:
                    self.ny = self.ny - 1
                    self.y = 0.5 * (y[1:] + y[:-1])
                    # left side of the cell
                    self.yl = y
                else:
                    self.y = y
                    self.yl = y
                if self.nz > 1:
                    self.nz = self.nz - 1
                    self.z = 0.5 * (z[1:] + z[:-1])
                    self.zl = z
                else:
                    self.z = z
                    self.zl = z
            elif point_type == "POINT_DATA":
                self.x = x
                self.y = y
                self.z = z

            grid_size = self.nx * self.ny * self.nz
            if grid_size != npoints:
                raise RuntimeError(
                    "Grid size ({}) is incompatible with number of points in the dataset ({})".format(
                        grid_size, npoints
                    )
                )
        else:
            # POLAR or SPHERICAL coordinates
            if self.nz == 1:
                is2d = 1
            else:
                is2d = 0

            s = fh.readline()  # POINTS NXNYNZ float
            slist = s.split()
            npoints = int(slist[1])
            points = np.fromfile(fh, dt, 3 * npoints)
            s = fh.readline()  # EXTRA LINE FEED

            self.points = points

            grid_size = self.nx * self.ny * self.nz
            if grid_size != npoints:
                raise RuntimeError(
                    "Grid size ({}) is incompatible with number of points in the dataset ({})".format(
                        grid_size, npoints
                    )
                )

            x1d = points[::3]
            y1d = points[1::3]
            z1d = points[2::3]

            xcart = np.transpose(x1d.reshape(self.nz, self.ny, self.nx))
            ycart = np.transpose(y1d.reshape(self.nz, self.ny, self.nx))
            zcart = np.transpose(z1d.reshape(self.nz, self.ny, self.nx))

            # Reconstruct the polar coordinate system
            if self.geometry == "polar":

                r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2)
                theta = np.unwrap(np.arctan2(ycart[0, :, 0], xcart[0, :, 0]))
                z = zcart[0, 0, :]

                s = fh.readline()  # CELL_DATA (NX-1)(NY-1)(NZ-1)
                slist = s.split()
                data_type = str(slist[0].decode("utf-8"))
                if data_type != "CELL_DATA":
                    print("ERROR: this routine expect CELL DATA as produced by idefix.")
                    fh.close()
                    return 0
                s = fh.readline()  # Line feed
                # Perform averaging on coordinate system to get cell centers
                # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
                if self.nx > 1:
                    self.nx = self.nx - 1
                    self.x = 0.5 * (r[1:] + r[:-1])
                    self.xl = r
                else:
                    self.x = r
                    self.xl = r
                if self.ny > 1:
                    self.ny = self.ny - 1
                    self.y = (0.5 * (theta[1:] + theta[:-1]) + np.pi) % (
                        2.0 * np.pi
                    ) - np.pi
                    self.yl = (theta + np.pi) % (2.0 * np.pi) - np.pi
                else:
                    self.y = theta
                    self.yl = theta
                if self.nz > 1:
                    self.nz = self.nz - 1
                    self.z = 0.5 * (z[1:] + z[:-1])
                    self.zl = z
                else:
                    self.z = z
                    self.zl = z

            # Reconstruct the spherical coordinate system
            if self.geometry == "spherical":
                if is2d:
                    r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2)
                    phi = np.unwrap(
                        np.arctan2(zcart[0, self.ny // 2, :], xcart[0, self.ny // 2, :])
                    )
                    theta = np.arccos(
                        ycart[0, :, 0]
                        / np.sqrt(xcart[0, :, 0] ** 2 + ycart[0, :, 0] ** 2)
                    )
                else:
                    r = np.sqrt(
                        xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2 + zcart[:, 0, 0] ** 2
                    )
                    phi = np.unwrap(
                        np.arctan2(
                            ycart[self.nx // 2, self.ny // 2, :],
                            xcart[self.nx // 2, self.ny // 2, :],
                        )
                    )
                    theta = np.arccos(
                        zcart[0, :, 0]
                        / np.sqrt(
                            xcart[0, :, 0] ** 2
                            + ycart[0, :, 0] ** 2
                            + zcart[0, :, 0] ** 2
                        )
                    )

                s = fh.readline()  # CELL_DATA (NX-1)(NY-1)(NZ-1)
                slist = s.split()
                data_type = str(slist[0].decode("utf-8"))
                if data_type != "CELL_DATA":
                    print("ERROR: this routine expect CELL DATA as produced by IDEFIX.")
                    fh.close()
                    return 0
                s = fh.readline()  # Line feed
                # Perform averaging on coordinate system to get cell centers
                # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
                if self.nx > 1:
                    self.nx = self.nx - 1
                    self.r = 0.5 * (r[1:] + r[:-1])
                    self.rl = r
                else:
                    self.r = r
                    self.rl = r
                if self.ny > 1:
                    self.ny = self.ny - 1
                    self.theta = 0.5 * (theta[1:] + theta[:-1])
                    self.thetal = theta
                else:
                    self.theta = theta
                    self.thetal = theta
                if self.nz > 1:
                    self.nz = self.nz - 1
                    self.phi = 0.5 * (phi[1:] + phi[:-1])
                    self.phil = phi
                else:
                    self.phi = phi
                    self.phil = phi

        ## From that point, the coordinates system is known.
        while 1:
            s = (
                fh.readline()
            )  # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
            # print repr(s)
            if len(s) < 2:  # leave if end of file
                break
            slist = s.split()
            datatype = str(slist[0].decode("utf-8"))
            varname = str(slist[1].decode("utf-8"))
            if datatype == "SCALARS":
                fh.readline()  # LOOKUP TABLE
                self.data[varname] = np.transpose(
                    np.fromfile(fh, dt, self.nx * self.ny * self.nz).reshape(
                        self.nz, self.ny, self.nx
                    )
                )
            elif datatype == "VECTORS":
                Q = np.fromfile(fh, dt, 3 * self.nx * self.ny * self.nz)

                self.data[varname + "_X"] = np.transpose(
                    Q[::3].reshape(self.nz, self.ny, self.nx)
                )
                self.data[varname + "_Y"] = np.transpose(
                    Q[1::3].reshape(self.nz, self.ny, self.nx)
                )
                self.data[varname + "_Z"] = np.transpose(
                    Q[2::3].reshape(self.nz, self.ny, self.nx)
                )

            else:
                raise RuntimeError("Unknown datatype '{}'".format(datatype))

            fh.readline()  # extra line feed

    def _load_particles(self, fh):
        raise NotImplementedError("Particles vtk are not supported yet !")

    def _setup_coordinates_from_native(self):
        if self.geometry == "spherical":
            native2attr = {
                "X1L_NATIVE_COORDINATES": "rl",
                "X1C_NATIVE_COORDINATES": "r",
                "X2L_NATIVE_COORDINATES": "thetal",
                "X2C_NATIVE_COORDINATES": "theta",
                "X3L_NATIVE_COORDINATES": "phil",
                "X3C_NATIVE_COORDINATES": "phi",
            }
        elif self.geometry in ("cartesian", "cylindrical", "polar"):
            native2attr = {
                "X1L_NATIVE_COORDINATES": "xl",
                "X1C_NATIVE_COORDINATES": "x",
                "X2L_NATIVE_COORDINATES": "yl",
                "X2C_NATIVE_COORDINATES": "y",
                "X3L_NATIVE_COORDINATES": "zl",
                "X3C_NATIVE_COORDINATES": "z",
            }

        for native_field, attr in native2attr.items():
            setattr(self, attr, self.native_coordinates[native_field])

    def __repr__(self):
        return "VTKDataset('%s')" % self.filename

# ////// public API //////
def readVTK(filename, geometry=None):
    r"""Read a VTK file for any geometry.

    This function is also compatible with particle VTK datasets

    With data older than Idefix 0.8, the *geometry* argument must be provided
    With more recent data it has no effect and can be ommited.
    """
    valid_geometries = tuple(KNOWN_GEOMETRIES.values())
    if geometry is not None and geometry not in valid_geometries:
        raise ValueError(
            "Unknown geometry '%s', expected one of %s" % (geometry, valid_geometries)
        )
    return VTKDataset(filename, geometry=geometry)


# Former geometry-specific readers (only for hydro datasets)
def readVTKCart(filename):
    warnings.warn(
        "the use of readVTKCart is discouraged. "
        "Use the generic readVTK function with geometry='cartesian'",
        stacklevel=2,
    )
    return readVTK(filename, geometry="cartesian")


# Read a vtk file
def readVTKPolar(filename):
    warnings.warn(
        "the use of readVTKPolar is discouraged. "
        "Use the generic readVTK function with geometry='polar'",
        stacklevel=2,
    )
    return readVTK(filename, geometry="polar")


# Read a vtk file
def readVTKSpherical(filename):
    warnings.warn(
        "the use of readVTKSpherical is discouraged. "
        "Use the generic readVTK function with geometry='spherical'",
        stacklevel=2,
    )
    return readVTK(filename, geometry="spherical")
