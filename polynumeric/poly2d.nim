import math

type
  TwoVarFunction* = proc (x, y: float): float

  Poly2d* = object
    coefs: seq[tuple[x, y: float]]

proc degree*(p: Poly2d): int =
  ## Returns the degree of the polynomial,
  ## that is the number of coefficients-1
  return p.coefs.len - 1

proc eval*(p: Poly2d, args: tuple[x, y: float]): float =
  ## Evaluates a polynomial function value for `x`
  ## quickly using Horners method
  var n = p.degree
  result = 0.0
  for i in countdown[int](n, 0):
    for j in countdown[int](n, 0):
      echo (i, j)
      let tmp = pow(args.x, i.toFloat) * p.coefs[i].x * pow(args.y, j.toFloat) * p.coefs[j].y
      result += tmp

proc eval*(p: Poly2d, x, y: float): float =
  p.eval((x: x, y: y))

proc `[]` *(p: Poly2d, idx: int): tuple[x, y: float] =
  ## Gets a coefficient of the polynomial.
  ## p[2] will returns the quadric term, p[3] the cubic etc.
  ## Out of bounds index will return 0.0.
  if idx < 0 or idx > p.degree:
    return (x: 0.0, y: 0.0)

  result.x = p.coefs[idx].x
  result.y = p.coefs[idx].y

proc `[]=` *(p: var Poly2d; idx: int, v: tuple[x, y: float]) =
  ## Sets an coefficient of the polynomial by index.
  ## p[2] set the quadric term, p[3] the cubic etc.
  ## If index is out of range for the coefficients,
  ## the polynomial grows to the smallest needed degree.
  assert(idx >= 0)

  if idx > p.degree: #polynomial must grow
    var oldlen = p.coefs.len
    p.coefs.setLen(idx + 1)
    for q in oldlen..<high(p.coefs):
      p.coefs[q].x = 0.0 # new-grown coefficients set to zero
      p.coefs[q].y = 0.0 # new-grown coefficients set to zero

  p.coefs[idx] = v

iterator items*(p: Poly2d): tuple[x, y: float] =
  ## Iterates through the coefficients of the polynomial.
  var i = p.degree
  while i >= 0:
    yield p[i]
    dec i


