# ported from https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/NormalToricVarieties/Divisors.m2
function monomials(R::MPolyRing, D::ToricDivisor)
  X = toric_variety(D)
  P = polyhedron(D)
  S = cox_ring(R, X)
  degs = matrix(ZZ, rays(X))
  coeffs = matrix(ZZ, n_rays(X), 1, coefficients(D))
  points = vec(lattice_points(P))
  # the stack blows up if we don't coerce the matrix element type from ZZRingElement to Int
  monomials = [S([1], [vec(Matrix{Int}(degs * p + coeffs))]) for p in points]
  return monomials
end

function monomials(D::ToricDivisor)
  X = toric_variety(D)
  P = polyhedron(D)
  S = cox_ring(X)
  degs = matrix(ZZ, rays(X))
  coeffs = matrix(ZZ, n_rays(X), 1, coefficients(D))
  points = vec(lattice_points(P))
  # the stack blows up if we don't coerce the matrix element type from ZZRingElement to Int
  monomials = [S([1], [vec(Matrix{Int}(degs * p + coeffs))]) for p in points]
  return monomials
end

function expected_dimension(s::Int, I::MPolyIdeal)
  n = ngens(base_ring(I)) - 1
  d = n - codim(I)
  return min(n, (s+1)*d + s)
end

function expected_dimension(s::Int, D::ToricDivisor)
  n = length(lattice_points(polyhedron(D))) - 1
  d = dim(variety(D))
  return min(n, (s+1)*d + s)
end

# Implementation based on definitiion of secant ideal via a change of coordinate on (s+1)-th
# embedding of X in P^n defined in Section 2.1 of "Equations of Secant Varieties: Geometry
# and Computation" by Jessica Sidman and Peter Vermeire.
# This reduces the number of variables in the ambient ring by one set of variables.
# Reference: https://doi.org/10.1007/978-3-642-19492-4_9
function secant_ideal(s::Int, I::MPolyIdeal)
  S = base_ring(I)
  K = coefficient_ring(S)
  n = ngens(S) - 1

  R, y = graded_polynomial_ring(K, "y" => (0:n, 1:s+1), cached=false)

  eta = hom(S, R, [y[j+1, s+1] - sum([y[j+1, i] for i in 1:s]) for j in 0:n])
  ring_homs = [hom(S, R, y[:, i]) for i in 1:s]
  J = sum(f -> f(I), ring_homs) + ideal(map(eta, gens(I)))
  psi = hom(R, S, vcat(zeros(S, s*(n+1)), gens(S)))

  return ideal(minimal_generating_set(psi(eliminate(J, vec(y[:, 1:s])))))
end

secant_ideal(s::Int, ntv::NormalToricVariety) = secant_ideal(s, toric_ideal(ntv))