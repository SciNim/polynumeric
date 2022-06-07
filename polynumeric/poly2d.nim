import arraymancer
import std/math

type
  TwoVarFunction* = proc (x, y: float): float

  Poly2d* = object
    coefs: Tensor[float]

proc initPoly2d*(coefs: Tensor[float]) : Poly2d =
  result.coefs = coefs

proc eval*(p: Poly2d, x, y: float): float =
  ## Evaluates a polynomial function value for `x`
  var
    nx = p.coefs.shape[0] - 1
    ny = p.coefs.shape[1] - 1
  result = 0.0
  for i in countdown[int](nx, 0):
    for j in countdown[int](ny, 0):
      result += pow(x, i.toFloat) * pow(y, j.toFloat) * (p.coefs[i, j])

proc eval*(p: Poly2d, args: tuple[x, y: float]): float =
  p.eval(args.x, args.y)
