import numpy as np
from scipy.interpolate import RegularGridInterpolator


def MakeNumpyFile():
    x = 2 ** np.arange(np.log2(1), np.log2(10), 0.2)
    y = 2 ** np.arange(np.log2(5), np.log2(10), 0.2)
    z = 2 ** np.arange(np.log2(2), np.log2(5), 0.2)

    print(x)
    print(y)
    print(z)
    xp, yp, zp = np.meshgrid(x, y, z, indexing="ij")

    data = xp + 2 * yp - zp

    np.save("x.npy", x)
    np.save("y.npy", y)
    np.save("z.npy", z)
    np.save("data.npy", data)
    f = RegularGridInterpolator((x, y, z), data)
    return f([2.7, 7.4, 3.9])
    # show the expected result
    #
    # print(f([2.7,7.4,3.9]))


if __name__ == "__main__":
    f = MakeNumpyFile()
    print("expected result:%f", f)
