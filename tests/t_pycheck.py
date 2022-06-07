import numpy as np
import numpy.polynomial.polynomial as poly

c = np.array([
        [-2.0, 3.0],
        [0.5, 1.2],
        [1.2, -1.3],
        [0.77, 0.88],
    ])

print(c)
res = poly.polyval2d(3.07, 4.05, c)
print(res)
