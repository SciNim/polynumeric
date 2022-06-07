## (c) Copyright 2013 Robert Persson

import std/math
import std/strutils

import ./polynumeric/poly2d

type
  OneVarFunction* = proc (x: float): float

  Poly* = object
    coefs: seq[float]

proc brent*(xmin, xmax: float, fn: OneVarFunction, tol: float, maxiter = 1000):
  tuple[rootx, rooty: float, success: bool] =
  ## Searches `fn` for a root between `xmin` and `xmax`
  ## using brents method. If the fn value at `xmin`and `xmax` has the
  ## same sign, `rootx`/`rooty` is set too the extrema value closest to x-axis
  ## and succes is set to false.
  ## Otherwise there exists at least one root and success is set to true.
  ## This root is searched for at most `maxiter` iterations.
  ## If `tol` tolerance is reached within `maxiter` iterations
  ## the root refinement stops and success=true.
  # see http://en.wikipedia.org/wiki/Brent%27s_method
  var
    a = xmin
    b = xmax
    c = a
    d = 1.0e308
    fa = fn(a)
    fb = fn(b)
    fc = fa
    s = 0.0
    fs = 0.0
    mflag: bool
    i = 0
    tmp2: float
  if fa*fb >= 0:
    if abs(fa) < abs(fb):
      return (a, fa, false)
    else:
      return (b, fb, false)
  if abs(fa) < abs(fb):
    swap(fa, fb)
    swap(a, b)
  while fb != 0.0 and abs(a - b) > tol:
    if fa != fc and fb != fc: # inverse quadratic interpolation
      s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb)
    else: #secant rule
      s = b - fb * (b - a) / (fb - fa)
    tmp2 = (3.0 * a + b) / 4.0
    if not((s > tmp2 and s < b) or (s < tmp2 and s > b)) or
      (mflag and abs(s - b) >= (abs(b - c) / 2.0)) or
      (not mflag and abs(s - b) >= abs(c - d) / 2.0):
      s = (a+b)/2.0
      mflag = true
    else:
      if (mflag and (abs(b - c) < tol)) or (not mflag and (abs(c - d) < tol)):
        s = (a+b) / 2.0
        mflag = true
      else:
        mflag = false
    fs = fn(s)
    d = c
    c = b
    fc = fb
    if fa * fs < 0.0:
      b = s
      fb = fs
    else:
      a = s
      fa = fs
    if abs(fa) < abs(fb):
      swap(a, b)
      swap(fa, fb)
    inc i
    if i > maxiter:
      break
  return (b, fb, true)

proc degree*(p: Poly): int =
  ## Returns the degree of the polynomial,
  ## that is the number of coefficients-1
  return p.coefs.len - 1

proc eval*(p: Poly, x: float): float =
  ## Evaluates a polynomial function value for `x`
  ## quickly using Horners method
  var n = p.degree
  result = p.coefs[n]
  dec n
  while n >= 0:
    result = result * x + p.coefs[n]
    dec n

proc `[]` *(p: Poly; idx: int): float =
  ## Gets a coefficient of the polynomial.
  ## p[2] will returns the quadric term, p[3] the cubic etc.
  ## Out of bounds index will return 0.0.
  if idx < 0 or idx > p.degree:
    return 0.0
  return p.coefs[idx]

proc `[]=` *(p: var Poly; idx: int, v: float) =
  ## Sets an coefficient of the polynomial by index.
  ## p[2] set the quadric term, p[3] the cubic etc.
  ## If index is out of range for the coefficients,
  ## the polynomial grows to the smallest needed degree.
  assert(idx >= 0)
  if idx > p.degree: #polynomial must grow
    var oldlen = p.coefs.len
    p.coefs.setLen(idx + 1)
    for q in oldlen..<high(p.coefs):
      p.coefs[q] = 0.0 #new-grown coefficients set to zero
  p.coefs[idx] = v

iterator items*(p: Poly): float =
  ## Iterates through the coefficients of the polynomial.
  var i = p.degree
  while i >= 0:
    yield p[i]
    dec i

proc clean*(p: var Poly; zerotol = 0.0) =
  ## Removes leading zero coefficients of the polynomial.
  ## An optional tolerance can be given for what's considered zero.
  var
    n = p.degree
    relen = false
  while n > 0 and abs(p[n]) <= zerotol: # >0 => keep at least one coefficient
    dec n
    relen = true
  if relen: p.coefs.setLen(n + 1)

proc `$`*(p: Poly): string =
  ## Gets a somewhat reasonable string representation of the polynomial
  ## The format should be compatible with most online function plotters,
  ## for example directly in google search
  result = ""
  var first = true #might skip + sign if first coefficient
  for idx in countdown(p.degree, 0):
    let a = p[idx]
    if a == 0.0:
      continue
    if a >= 0.0 and not first:
      result.add('+')
    first = false
    if a != 1.0 or idx == 0:
      result.add(formatFloat(a, ffDefault, 0))
    if idx >= 2:
      result.add("x^" & $idx)
    elif idx == 1:
      result.add("x")
  if result == "":
    result = "0"

proc derivative*(p: Poly): Poly =
  ## Returns a new polynomial, which is the derivative of `p`
  newSeq[float](result.coefs, p.degree)
  for idx in 0..high(result.coefs):
    result.coefs[idx] = p.coefs[idx + 1] * float(idx + 1)

proc diff*(p: Poly, x: float): float =
  ## Evaluates the differentiation of a polynomial with
  ## respect to `x` quickly using a modifed Horners method
  var n = p.degree
  result = p[n] * float(n)
  dec n
  while n >= 1:
    result = result * x + p[n] * float(n)
    dec n

proc integral*(p: Poly): Poly =
  ## Returns a new polynomial which is the indefinite
  ## integral of `p`. The constant term is set to 0.0
  newSeq(result.coefs, p.coefs.len + 1)
  result.coefs[0] = 0.0 #constant arbitrary term, use 0.0
  for i in 1..high(result.coefs):
    result.coefs[i] = p.coefs[i-1] / float(i)

proc integrate*(p: Poly; xmin, xmax: float): float =
  ## Computes the definite integral of `p` between `xmin` and `xmax`
  ## quickly using a modified version of Horners method
  var
    n = p.degree
    s1 = p[n] / float(n + 1)
    s2 = s1
    fac: float
  dec n
  while n >= 0:
    fac = p[n] / float(n + 1)
    s1 = s1 * xmin + fac
    s2 = s2 * xmax + fac
    dec n
  result = s2 * xmax - s1 * xmin

proc initPoly*(cofs: varargs[float]): Poly =
  ## Initializes a polynomial with given coefficients.
  ## The most significant coefficient is first, so to create x^2-2x+3:
  ## intiPoly(1.0,-2.0,3.0)
  if len(cofs) <= 0:
    result.coefs = @[0.0] #need at least one coefficient
  else:
    # reverse order of coefficients so indexing matches degree of
    # coefficient...
    result.coefs = @[]
    for idx in countdown(cofs.len-1, 0):
      result.coefs.add(cofs[idx])
  result.clean #remove leading zero terms

proc divMod*(p, d: Poly; q, r: var Poly) =
  ## Divides `p` with `d`, and stores the quotinent in `q` and
  ## the remainder in `d`
  var
    pdeg = p.degree
    ddeg = d.degree
    power = p.degree - d.degree
    ratio: float
  r.coefs = p.coefs #initial remainder=numerator
  if power < 0: #denominator is larger than numerator
    q.coefs = @[0.0] #quotinent is 0.0
    return # keep remainder as numerator
  q.coefs = newSeq[float](power + 1)
  for i in countdown(pdeg, ddeg):
    ratio = r.coefs[i] / d.coefs[ddeg]
    q.coefs[i-ddeg] = ratio
    r.coefs[i] = 0.0
    for j in 0 ..< ddeg:
      var idx = i-ddeg+j
      r.coefs[idx] = r.coefs[idx] - d.coefs[j] * ratio
  r.clean # drop zero coefficients in remainder

proc `+` *(p1: Poly, p2: Poly): Poly =
  ## Adds two polynomials
  var n = max(p1.coefs.len, p2.coefs.len)
  newSeq(result.coefs, n)
  for idx in countup(0, n - 1):
    result[idx] = p1[idx] + p2[idx]
  result.clean # drop zero coefficients in remainder

proc `*` *(p1, p2: Poly): Poly =
  ## Multiplies the polynomial `p1` with `p2`
  var
    d1 = p1.degree
    d2 = p2.degree
    n = d1 + d2
    idx: int
  newSeq(result.coefs, n)
  for i1 in countup(0, 1):
    for i2 in countup(0, d2):
      idx = i1 + i2
      result[idx] = result[idx] + p1[i1] * p2[i2]
  result.clean

proc `*` *(p: Poly, f: float): Poly =
  ## Multiplies the polynomial `p` with a real number
  newSeq(result.coefs, p.coefs.len)
  for i in 0..high(p.coefs):
    result[i] = p.coefs[i] * f
  result.clean

proc `*` *(f: float, p: Poly): Poly =
  ## Multiplies a real number with a polynomial
  result = p * f

proc `-`*(p: Poly): Poly =
  ## Negates a polynomial
  result = p
  for i in 0 ..< result.coefs.len:
    result.coefs[i] = -result.coefs[i]

proc `-`*(p1, p2: Poly): Poly =
  ## Subtract `p1` with `p2`
  var n = max(p1.coefs.len, p2.coefs.len)
  newSeq(result.coefs, n)
  for idx in countup(0, n-1):
    result[idx] = p1[idx] - p2[idx]
  result.clean # drop zero coefficients in remainder

proc `/`*(p: Poly, f: float): Poly =
  ## Divides polynomial `p` with a real number `f`
  newSeq(result.coefs, p.coefs.len)
  for i in 0..high(p.coefs):
    result[i] = p.coefs[i] / f
  result.clean

proc `/` *(p, q: Poly): Poly =
  ## Divides polynomial `p` with polynomial `q`
  var dummy: Poly
  p.divMod(q, result, dummy)

proc `mod` *(p, q: Poly): Poly =
  ## Computes the polynomial modulo operation,
  ## that is the remainder of `p`/`q`
  var dummy: Poly
  p.divMod(q, dummy, result)

proc normalize*(p: var Poly) =
  ## Multiplies the polynomial inplace by a term so that
  ## the leading term is 1.0.
  ## This might lead to an unstable polynomial
  ## if the leading term is zero.
  p = p / p[p.degree]

proc solveQuadric*(a, b, c: float; zerotol = 0.0): seq[float] =
  ## Solves the quadric equation `ax^2+bx+c`, with a possible
  ## tolerance `zerotol` to find roots of curves just 'touching'
  ## the x axis. Returns sequence with 0,1 or 2 solutions.
  var p, q, d: float
  p = b / (2.0 * a)
  if p == Inf or p == NegInf: #linear equation..
    var linrt = -c/b
    if linrt == Inf or linrt == NegInf: #constant only
      return @[]
    return @[linrt]
  q = c / a
  d = p*p-q
  if d < 0.0:
    #check for inside zerotol range for neg. roots
    var err = a * p * p - b * p + c #evaluate error at parabola center axis
    if err <= zerotol: return @[-p]
    return @[]
  else:
    var sr = sqrt(d)
    result = @[-sr-p, sr-p]

proc getRangeForRoots(p: Poly): tuple[xmin, xmax: float] =
  ## helper function for `roots` function
  ## quickly computes a range, guaranteed to contain
  ## all the real roots of the polynomial
  # see http://www.mathsisfun.com/algebra/polynomials-bounds-zeros.html
  var
    deg = p.degree
    d = p[deg]
    bound1, bound2: float
  for i in countup(0, deg):
    var c = abs(p.coefs[i] / d)
    bound1 = max(bound1, c+1.0)
    bound2 = bound2+c
  bound2 = max(1.0, bound2)
  result.xmax = min(bound1, bound2)
  result.xmin = -result.xmax

proc addRoot(p: Poly, res: var seq[float], xp0, xp1, tol, zerotol, mergetol: float, maxiter: int) =
  ## helper function for `roots` function. Try to do a numeric search for a single root
  ## in range xp0-xp1, adding it to `res` (allocating `res` if nil)
  var br = brent(xp0, xp1, proc(x: float): float = p.eval(x), tol)
  if br.success:
    if res.len == 0 or br.rootx >= res[high(res)] + mergetol: # don' t add equal roots.
      res.add(br.rootx)
  else:
    #this might be a 'touching' case, check function value against zero tolerance
    if abs(br.rooty) <= zerotol:
      if res.len == 0 or br.rootx >= res[high(res)] + mergetol: # don't add equal roots.
        res.add(br.rootx)

proc roots*(p: Poly, tol = 1.0e-9, zerotol = 1.0e-6, mergetol = 1.0e-12, maxiter = 1000): seq[float] =
  ## Computes the real roots of the polynomial `p`
  ## `tol` is the tolerance used to break searching for each root when reached.
  ## `zerotol` is the tolerance, which is 'close enough' to zero to be considered a root
  ## and is used to find roots for curves that only 'touch' the x-axis.
  ## `mergetol` is the tolerance, of which two x-values are considered being the same root.
  ## `maxiter` can be used to limit the number of iterations for each root.
  ## Returns a (possibly empty) sorted sequence with the solutions.
  var deg = p.degree
  if deg <= 0: #constant only => no roots
    return @[]
  elif p.degree == 1: #linear
    var linrt = -p.coefs[0] / p.coefs[1]
    if linrt == Inf or linrt == NegInf:
      return @[] #constant only => no roots
    return @[linrt]
  elif p.degree == 2:
    return solveQuadric(p.coefs[2], p.coefs[1], p.coefs[0], zerotol)
  else:
    # degree >=3 , find min/max points of polynomial with recursive
    # derivative and do a numerical search for root between each min/max
    var rng = p.getRangeForRoots()
    var minmax = p.derivative.roots(tol, zerotol, mergetol)
    result = @[]
    if minmax.len > 0: #ie. we have minimas/maximas in this function
      for x in minmax.items:
        addRoot(p, result, rng.xmin, x, tol, zerotol, mergetol, maxiter)
        rng.xmin = x
    addRoot(p, result, rng.xmin, rng.xmax, tol, zerotol, mergetol, maxiter)

# imported here as it's ``only`` used for the fit
import arraymancer / [tensor, linear_algebra]
proc polyFit*[T: seq[float] | Tensor[float]](x, y: T, polyOrder: int): Tensor[float] =
  ## Performs a linear least squares fit to the data x and y with a polynomial
  ## of order `polyOrder`.
  ##
  ## NOTE: As it uses LAPACK's least squares solver, it does depend on LAPACK!
  ## See the arraymancer README for the correct compilation flags for your system.
  when T is seq[float]:
    let x = x.toTensor
    let y = y.toTensor

  # actual poly order is + 1, because `0` corresponds to polynomial of order 0, a constant
  var A = vandermonde(x, polyOrder + 1)
  # scale lhs to improve condition number and solve
  let scale = sqrt((A *. A).sum(axis = 0))
  A = A /. scale
  var (coeffs, resids, rank, s) = least_squares_solver(A, y)
  coeffs = (coeffs /. scale.squeeze)  # broadcast scale coefficients
  result = coeffs
