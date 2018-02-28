import numpy as np

def compute_lattice_constants(bsize, tilt_factors):

    for key in ["x", "y", "z"]:
        if key.lower() not in bsize and key.upper() not in bsize:
            raise KeyError("Could not find key '%s'." % key)

    for key in ["xy", "xz", "yz"]:
        if key.lower() not in tilt_factors and key.upper() not in tilt_factors:
            raise KeyError("Could not find key '%s'." % key)

    lx = bsize['x']
    ly = bsize['y']
    lz = bsize['z']

    xy = tilt_factors['xy']
    xz = tilt_factors['xz']
    yz = tilt_factors['yz']

    a = lx
    b = np.sqrt(np.power(ly, 2) + np.power(xy, 2))
    c = np.sqrt(np.power(lz, 2) + np.power(xz, 2) + np.power(yz, 2))

    if np.isclose(b * c, 0.0):
        raise ZeroDivisionError("One of the box sizes is zero")

    cos_alpha = (xy *  xz + ly * yz) / (b * c)
    cos_beta = xz / c
    cos_gamma = xy / b

    alpha = np.arccos(cos_alpha)
    beta = np.arccos(cos_beta)
    gamma = np.arccos(cos_gamma)

    return {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

