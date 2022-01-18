import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class Anisotropic:
    def __init__(self, a, C):
        a = np.asarray(a, dtype=np.float64)
        C = np.asarray(C, dtype=np.float64)
        assert a.shape == (3,)
        assert C.shape == (3, 3)

        # A = a a^T 
        self.A = (a).reshape(3, 1) @ a.T
        # `$I_5$` = tr(CA)
        self.I_5 = np.trace(C @ self.A)
        # `$\frac{∂²I₅}{∂f²}$` = 2[A₁,₁I₃  A₁,₂I₃  A₁,₃I₃
    #                A₂,₁I₃  A₂,₂I₃  A₂,₃I₃
    #                A₃,₁I₃  A₃,₂I₃  A₃,₃I₃] 
        frac_partial_differential_2I5_partial_differential_f2_0 = np.block([[self.A[1-1, 1-1] * np.identity(3), self.A[1-1, 2-1] * np.identity(3), self.A[1-1, 3-1] * np.identity(3)], [self.A[2-1, 1-1] * np.identity(3), self.A[2-1, 2-1] * np.identity(3), self.A[2-1, 3-1] * np.identity(3)], [self.A[3-1, 1-1] * np.identity(3), self.A[3-1, 2-1] * np.identity(3), self.A[3-1, 3-1] * np.identity(3)]])
        self.frac_partial_differential_2I5_partial_differential_f2 = 2 * frac_partial_differential_2I5_partial_differential_f2_0

