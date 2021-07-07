import unittest

import arraymancer
import ../polynumeric/poly2d

proc test() =
  var coefs = [[-2.0, 3.0],
          [0.5, 1.2],
          [1.2, -1.3],
          [0.77, 0.88],
          ].toTensor

  var p = initPoly2d(coefs)
  block:
    var res = p.eval((x: 3.07, y: 4.05))
    echo res

  block:
    var res = p.eval(3.07, 4.05)
    echo res

when isMainModule:
  test()
