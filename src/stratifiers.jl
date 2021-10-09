# stratifications.jl
# Paper-specific stratification functions.

"""
US Census Region
"""
function regionate!(cc::AbstractCICModel; regiontype = :region, datpath = false)
  if !(regiontype .∈ Ref([:region, :state, :division]))
    error("supply valid region type")
  end

  cc.stratifier = :Region;

  if !datpath
    creg = CSV.read(download("https://raw.githubusercontent.com/kjhealy/fips-codes/master/county_fips_master.csv"), DataFrame)
    @subset!(creg, $regiontype .!= "NA")
  else
    creg = read(datpath, DataFrame)
    @subset!(creg, $regiontype .!= "NA")
  end

  rnme = Symbol(string(regiontype) * "_name");
  X = sort(unique(creg[!, [regiontype, rnme]]), regiontype)[!, rnme];

  dd = Dict{Int, Int}();
  sizehint!(dd, nrow(creg))
  for (i, v) in enumerate(creg[!, regiontype])
    vx = tryparse(Int, v)
    dd[creg[i, :fips]] = isnothing(vx) ? 0 : vx
  end

  cc.matches[!, :stratum] = [dd[cc.matches[i, :treatunit]] for i in 1:nrow(cc.matches)];

  cc.meanbalances[!, :stratum] = [dd[cc.balances[i, :treatunit]] for i in 1:nrow(cc.meanbalances)];
    
  return cc, Dict(1:length(X) .=> X)
end

"""
Treatment Date
"""
function datestrat!(
  cc::AbstractCICModel;
  qtes = [0, 0.25, 0.5, 0.75, 1.0],
  datestart = Date("2020-03-01")
)

  cc.stratifier = Symbol("Primary Date");

  tts = unique(cc.matches[!, :treattime]);

  udtes = unique(tts); # we want quantiles of unique
  # udtes = unique(@view dat[dat[!, treatment] .== 1, :date]);

  X = sort(Int.(round.(quantile(udtes, qtes), digits = 0)));
  Xlen = length(X);

  cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
  @eachrow! cc.matches begin
    :stratum = assignq(:treattime, X, Xlen)
  end

  cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));
  @eachrow! cc.meanbalances begin
    :stratum = assignq(:treattime, X, Xlen)
  end

  dteq = string.(datestart + Day.(X));

  return cc, label_variablestrat(dteq)
end

# Population Density
# cc, labs = variablestrat!(
#   cc, dat, Symbol("Pop. Density");
#   qtes = [0, 0.25, 0.5, 0.75, 1.0]
# )

# Treatment Date - Date of First Case
function primarydistancestrat!(cc, dat; qtes = [0, 0.25, 0.5, 0.75, 1.0])

  cc.stratifier = Symbol("Date of First Case to Primary");

  # get first case dates
  uid = unique(dat[!, :fips]);
  
  dd = Dict{Int64, Int64}();
  sizehint!(dd, length(uid));
  # mm = maximum(dat[!, :running])
  for u in uid
    udat = @view dat[dat[!, :fips] .== u, [:fips, :running, :casescum]];
    udat = sort(udat, [:fips, :running], view = true);
    ff = findfirst(udat[!, :casescum] .> 0);
    if !isnothing(ff)
      dd[u] = ff
    end
  end

  # X = sort(quantile(values(dd), qtes));
  X = sort(Int.(round.(quantile(values(dd), qtes), digits = 0)))
  Xlen = length(X);

  cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
  @eachrow! cc.matches begin
    :stratum = assignq(udict[:treatunit], X, Xlen)
  end

  cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));
  @eachrow! cc.balances begin
    :stratum = assignq(udict[:treatunit], X, Xlen)
  end

  return cc, label_variablestrat(string.(X))
end

# Trump's Share of the Vote in 2016
# variablestrat!(
#   cc, dat, Symbol("Trump 2016 Vote Share");
#   qtes = [0, 0.25, 0.5, 0.75, 1.0]
# )

# # Voter Turnout
# variablestrat!(
#   cc, dat, Symbol("In-person Turnout Rate");
#   qtes = [0, 0.25, 0.5, 0.75, 1.0]
# )

# # Cumulative Case Rate
# variablestrat!(
#   cc, dat, Symbol("Cum. Case Rate");
#   qtes = [0, 0.25, 0.5, 0.75, 1.0], zerosep = false
# )

# # Cumulative Death Rate
# variablestrat!(
#   cc, dat, Symbol("Cum. Death Rate");
#   qtes = [0, 0.25, 0.5, 0.75, 1.0], zerosep = false
# )
