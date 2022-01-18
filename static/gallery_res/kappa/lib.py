import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class kappa:
    def __init__(self, a, t, θ):
        assert np.ndim(a) == 0
        assert np.ndim(t) == 0
        assert np.ndim(θ) == 0

        self.θ = θ
        self.a = a
        # x = a cos(t) cot(t)
        self.x = a * np.cos(t) * 1/np.tan(t)
        # y = a cos(t)
        self.y = a * np.cos(t)
        # r = atan(θ)
        self.r = np.arctan(θ)

    def κ(self, θ):
        assert np.ndim(θ) == 0

        return (8 * (3 - np.power(np.sin(θ), 2)) * np.power(np.sin(θ), 4)) / (a * np.power((np.power(np.sin(2 * θ), 2) + 4), (3 / 2)))

    def phi(self, θ):
        assert np.ndim(θ) == 0

        return -np.arctan(1 / 2 * np.sin(2 * θ))

