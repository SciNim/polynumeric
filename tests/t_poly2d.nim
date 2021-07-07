import unittest
import ../polynumeric/poly2d


proc test() =
  var p_coefs = @[
    (x: -2.0, y: 3.0),
    (x: 0.5, y: 1.2),
    (x: 1.2, y: -1.3),
  ]

  var p : Poly2d
  p[0] = (x: -2.0, y: 3.0)
  p[1] = (x: 0.5, y: 1.2)
  p[2] = (x: 1.2, y: -1.3)
  check p[0] == (x: -2.0, y: 3.0)
  check p[1] == (x: 0.5, y: 1.2)
  check p[2] == (x: 1.2, y: -1.3)

  p[3] = (x: 0.77, y: 0.88)
  check p[3] == (x: 0.77, y: 0.88)

  block:
    var res = p.eval((x: 3.07, y: 4.05))
    echo res

  block:
    var res = p.eval(3.07, 4.05)
    echo res

when isMainModule:
  test()
