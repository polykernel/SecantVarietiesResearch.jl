# ported from https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/NormalToricVarieties/Divisors.m2
function monomials(R::MPolyRing, D::ToricDivisor)
  X = toric_variety(D)
  P = polyhedron(D)
  S = cox_ring(R, X)
  degs = matrix(ZZ, rays(X))
  coeffs = matrix(coefficients(D))
  lat_points = lattice_points(P)
  # the stack blows up if we don't coerce the matrix element type from ZZRingElement to Int
  monomials = [S([1], [vec(collect(Int, degs * p + coeffs))]) for p in lat_points]
  return monomials
end

function monomials(D::ToricDivisor)
  X = toric_variety(D)
  S, _ = polynomial_ring(coefficient_ring(X), coordinate_names(X), cached=false)
  return monomials(S, D)
end

function expected_dimension(s::Int, I::MPolyIdeal)
  n = ngens(base_ring(I)) - 1
  d = n - codim(I)
  return min(n, (s+1)*d + s)
end

function expected_dimension(s::Int, D::ToricDivisor)
  n = length(lattice_points(polyhedron(D))) - 1
  d = dim(toric_variety(D))
  return min(n, (s+1)*d + s)
end

function vanishing_ideal(K::Ring, D::ToricDivisor)
  X = toric_variety(D)

  R, r_gens = graded_polynomial_ring(K, ["x_$i" for i in 0:length(coordinate_names(X))-1], cached=false)

  L = monomials(R, D)

  BR = parent(L[1])
  mo = neginvlex(BR)
  sort!(L, lt=(x, y) -> cmp(mo, x, y) == -1)

  S, _ = graded_polynomial_ring(K, ["x_$i" for i in 0:length(L)-1], cached=false)

  phi = hom(S, R, [evaluate(l, r_gens) for l in L])
  I = kernel(phi)

  return I
end

function vanishing_ideal(D::ToricDivisor)
  X = toric_variety(D)
  R = coefficient_ring(X)
  return vanishing_ideal(R, D)
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

secant_ideal(K::Ring, s::Int, D::ToricDivisor) = secant_ideal(s, vanishing_ideal(K, D))
secant_ideal(s::Int, D::ToricDivisor) = secant_ideal(s, vanishing_ideal(D))


_random_fieldelem(KK::RealField, args...) = rand(KK, args...)
_random_fieldelem(KK::T, args...) where T <: Union{FqField, FracField} = rand(make(KK, args...))

function terracini_dimension(s::Int, D::ToricDivisor, KK::Union{FqField, FracField, RealField}, args...)
  d = dim(toric_variety(D))
  P = polyhedron(D)
  points = Int.(reduce(hcat, lattice_points(P)))
  A = begin
    c = ncols(points)
    min_point = [minimum(point) for point in eachrow(points)]
    T = reduce(hcat, fill(min_point, c))
    mat = points - T
    col_sum = sum.(eachcol(mat))
    m = maximum(col_sum)
    vcat(reshape(fill(m, c) - col_sum, 1, :), mat)
  end
  R, _ = graded_polynomial_ring(KK, ["x_$i" for i in 0:d], cached=false)
  M = [R([1], [copy(exp)]) for exp in eachcol(A)]
  random_points = [[_random_fieldelem(KK, args...) for _ in 0:d] for _ in 0:s]
  jac_mat = jacobian_matrix(M)
  tangent_span = reduce(vcat, [map_entries(p -> evaluate(p, points), jac_mat) for points in random_points])
  return rank(tangent_span) - 1
end
