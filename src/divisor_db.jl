const SMALL_SMOOTH_PROJETIVE_DIVISORS_FILE = joinpath(
  DATA_DIR, "small_smooth_projective_toric_varieties.ndjson"
)

const SMOOTH_FANO_TORIC_VARIETIES_FILES = [
  joinpath(DATA_DIR, "smooth_fano_toric_varieties.ndjson"),
  joinpath(DATA_DIR, "smooth_fano_toric_varieties_5.ndjson"),
  joinpath(DATA_DIR, "smooth_fano_toric_varieties_6.ndjson"),
]

const SMALL_SMOOTH_PROJETIVE_DIVISORS_DB = begin
  hashtable = Dict()

  open(SMALL_SMOOTH_PROJETIVE_DIVISORS_FILE) do f
    for ln in eachline(f)
      data = JSON3.read(ln)
      hashtable[(data[1], data[2])] = Dict(
        :rays => copy(data[3]),
        :max_cones => broadcast(x -> 1 .+ x, copy(data[4])),
        :coefficients => copy(data[5]),
      )
    end
  end

  hashtable
end

const SMOOTH_FANO_TORIC_VARIETIES_DB = begin
  hashtable = Dict()

  for file in SMOOTH_FANO_TORIC_VARIETIES_FILES
    open(file) do f
      for ln in eachline(f)
        data = JSON3.read(ln)
        hashtable[(data[1], data[2])] = Dict(
          :rays => copy(data[3]), :max_cones => broadcast(x -> 1 .+ x, copy(data[4]))
        )
      end
    end
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

# ported from https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/NormalToricVarieties/Divisors.m2
function smooth_fano_toric_variety(d::Int, i::Int)
  key = (d, i)
  if !haskey(SMOOTH_FANO_TORIC_VARIETIES_DB, key)
    error("Index $key is not in the database")
  end
  variety_data = SMOOTH_FANO_TORIC_VARIETIES_DB[(d, i)]
  # the non_redundant flag is needed to preserve the order of the rays
  X = normal_toric_variety(
    IncidenceMatrix(variety_data[:max_cones]), variety_data[:rays]; non_redundant=true
  )
  return X
end
