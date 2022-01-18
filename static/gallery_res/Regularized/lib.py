import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class Regularized:
    def __init__(self, f, s, F, r_ε, a, b, ε):
        f = np.asarray(f, dtype=np.float64)
        F = np.asarray(F, dtype=np.float64)
        assert f.shape == (3,)
        assert np.ndim(s) == 0
        assert F.shape == (3, 3)
        assert np.ndim(r_ε) == 0
        assert np.ndim(a) == 0
        assert np.ndim(b) == 0
        assert np.ndim(ε) == 0

        self.r_ε = r_ε
        self.a = a
        self.b = b
        self.ε = ε
        self.f = f
        self.F = F
        self.s = s
        # placeholder = 1
        self.placeholder = 1

    def row_ε(self, r):
        r = np.asarray(r, dtype=np.float64)
        assert r.shape == (3,)

        return (15 * r_ε / 8 + 1 / np.power(r_ε, 3))

    def u_ε(self, r):
        r = np.asarray(r, dtype=np.float64)
        assert r.shape == (3,)

        return ((a - b) / r_ε * np.identity(3) + (b / np.power(r_ε, 3) * r).reshape(3, 1) @ r.T + a / 2 * np.power(ε, 2) / np.power(r_ε, 3) * np.identity(3)) @ f

    def t_ε(self, r):
        r = np.asarray(r, dtype=np.float64)
        assert r.shape == (3,)

        return -a * (1 / np.power(r_ε, 3) + 3 * np.power(ε, 2) / (2 * np.power(r_ε, 5))) * F @ r

    def s_ε(self, r):
        r = np.asarray(r, dtype=np.float64)
        assert r.shape == (3,)

        return (2 * b - a) * (1 / np.power(r_ε, 3) + 3 * np.power(ε, 2) / (2 * np.power(r_ε, 5))) * (s * r)

    def p_ε(self, r):
        r = np.asarray(r, dtype=np.float64)
        assert r.shape == (3,)

        return (2 * b - a) / np.power(r_ε, 3) * F @ r - 3 / (2 * np.power(r_ε, 5)) * (2 * b * ((r.T.reshape(1, 3) @ F @ r).item()) * np.identity(3) + a * np.power(ε, 2) * F) @ r

