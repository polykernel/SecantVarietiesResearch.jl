const SMALL_SMOOTH_PROJETIVE_DIVISORS_FILE = joinpath(
  DATA_DIR, "small_smooth_projective_toric_varieties.ndjson"
)

const SMALL_SMOOTH_PROJETIVE_DIVISORS_DB = open(SMALL_SMOOTH_PROJETIVE_DIVISORS_FILE) do f
  hashtable = Dict()
  for ln in eachline(f)
    data = JSON3.read(ln)
    hashtable[(data[1], data[2])] = Dict(
      :rays => copy(data[3]),
      :max_cones => broadcast(x -> 1 .+ x, copy(data[4])), # julia arrays/matrices are 1-indexed
      :coefficients => copy(data[5]),
    )
  end
  hashtable
end

# ported from https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/NormalToricVarieties/Divisors.m2
function small_ample_toric_divisor(d::Int, i::Int)
  key = (d, i)
  if !haskey(SMALL_SMOOTH_PROJETIVE_DIVISORS_DB, key)
    error("Index $key is not in the database")
  end
  divisor_data = SMALL_SMOOTH_PROJETIVE_DIVISORS_DB[(d, i)]
  # the non_redundant flag is needed to preserve the order of the rays
  X = normal_toric_variety(
    IncidenceMatrix(divisor_data[:max_cones]), divisor_data[:rays]; non_redundant=true
  )
  D = toric_divisor(X, divisor_data[:coefficients])
  return D
end