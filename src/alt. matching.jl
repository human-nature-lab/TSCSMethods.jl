# new new matching
# nested DataFrame


"""
inverted covariance matrix for mahalanobis distance (all units at t)

inverted sqrt(vars) for balance score calculations (treated units at t)
"""
function calculate_sample_Σs!(
  ut, Σinvdict, σdict, dt, c_treatment, cmat,
  variancesonly::Bool
)
  
  for i = eachindex(ut)
    uti = ut[i]
    c_t = dt .== uti
    Σ = cov(cmat[c_t, :])

    if variancesonly
      Σ[Not(diagind(Σ))] .= 0
    end

    c2 = c_t .& c_treatment;
    if (sum(c2) > 1)
      Σ_treated = cov(cmat[c2, :])
      σdict[uti] = 1 ./ sqrt.(diag(Σ_treated)) # for balance score
    end
    Σinvdict[uti] = pinv(Σ)
  end
  return Σinvdict, σdict
end

function get_Σs(
  dat, uid, covariates, t, id, treatment;
  variancesonly = true
)

  sdat = dat[!, vcat(t, covariates, id, treatment)]

  ut = unique(sdat[:, t])
  cmat = Matrix(sdat[!, covariates])

  c_treatment = sdat[!, id] .∈ Ref(uid)

  Σinvdict = Dict{Int64, Matrix{Float64}}();
  σdict = Dict{Int64, Vector{Float64}}();

  calculate_sample_Σs!(
    ut, Σinvdict, σdict, sdat[!, t], c_treatment, cmat,
    variancesonly
  )
  return Σinvdict, σdict
end

function distance_sums!(
  tmin, Σinvdict, covariates, c1tusdf, c1ousdf, LF, Lrnge, caldists, mdist, lcnt
)
  
  for l = eachindex(Lrnge) # mmin:mmax
    ltime = LF[l]
    if ltime >= tmin
      lcnt[1] += 1

      il = Vector{Float64}(undef, length(covariates))
      jl = similar(il)
      
      # fillrow!(il, jl, covariates, c1tusdf, c1ousdf, lcnt);

      fillrow!(
        il, jl, covariates, c1tusdf, c1ousdf, lcnt
      )

      mdist[1] += mahalanobis(il, jl, Σinvdict[LF[l]])
      distance_sums_inner!(caldists, covariates, il, jl, Σinvdict, LF, l)
      
    end
  end
  return lcnt, mdist, caldists
end

function fillrow!(il, jl, covariates, c1tusdf, c1ousdf, lcnt)
  for (cnum, covar) = enumerate(covariates)
    il[cnum] = c1tusdf[lcnt[1], covar] # row for iu at l
    jl[cnum] = c1ousdf[lcnt[1], covar] # row for ju at l

  end
return il,jl
end

function distance_sums_inner!(caldists, covariates, il, jl, Σinvdict, LF, l)
    for cnum = eachindex(covariates)
        caldists[cnum] += weuclidean(
          il[cnum], jl[cnum],
          Σinvdict[LF[l]][cnum,cnum]
        )
    end
    return caldists
end

function setupmatching(cc::cicmodel, dat::DataFrame)

  
  tmin = minimum(dat[!, t]);

  ttobs = unique(dat[dat[!, cc.treatment] .== 1, [cc.t, cc.id]]);

  import IterTools
  
  CP = IterTools.product(
    unique(dat[!, cc.id]),
    [(ttobs.running[i], ttobs[i, id]) for i in 1:nrow(ttobs)]
  );

  cc.matches = DataFrame(
    ttime = zeros(Int64, length(CP)),
    treatunit = zeros(Int64, length(CP)),
    matchunit = zeros(Int64, length(CP))
  );

  treatpretreat = Vector{Float64}(undef, size(collect(CP), 2));

  matchpretreat = Matrix{Float64}(undef, size(collect(CP), 1), size(collect(CP), 2));

  treatmatchperavg = Matrix{Vector{Union{Float64, Missing}}}(undef, size(collect(CP), 1), size(collect(CP), 2));

  matchmatchperavg = similar(treatmatchperavg);
  treatoutcome = similar(treatmatchperavg);
  matchoutcome = similar(treatmatchperavg);
  mdistes = similar(treatmatchperavg);
  possible = Matrix{Vector{Bool}}(undef, size(collect(CP), 1), size(collect(CP), 2));

  cals = [similar(treatmatchperavg) for cv in 1:length(cc.covariates)];

  #=
  1158.105712 seconds (3.90 G allocations: 821.754 GiB, 52.67% gc time, 0.00% compilation time)
  
  ~20 minutes
  =#

  lK = length(cc.fmin:cc.fmax)
  lM = length(mmin:mmax);

  @time @inbounds Threads.@threads for c = eachindex(1:size(CP)[2])

    cp = collect(CP)[:,  c];
    tu = cp[1][2][2]
    tt = cp[1][2][1]

    ttl = tt + cc.fmin + cc.mmin;
    ttu = tt + cc.fmax;
    sdf = @view dat[(dat[!, cc.t] .>= ttl) .& (dat[!, cc.t] .<= ttu), :];
    tusdf = @view sdf[sdf[!, cc.id] .== tu, :];

    LF = (tt + cc.fmin) + cc.mmin : (tt + cc.fmax) + cc.mmax;
  
    tstrt = max((tt + cc.fmin) + cc.mmin, tmin) # begin of data or earliest L
    # this is the first time point in sdf
    # assumes same start date for all units
    trtdx = tt - tstrt + 1; # first index is 1

    if trtdx < 2 # treatment must be present, and day before treatment
      # possibles[ijk, :poss] = false; # default already
      continue
    end

    treatpretreat[c] = tusdf[trtdx + cc.reference, cc.outcome];
    
    for r = eachindex(1:size(CP)[1])
      mu = cp[r][1]
      judat = @view sdf[sdf[!, cc.id] .== mu, :];

      matchpretreat[r, c] = judat[trtdx + cc.reference, cc.outcome];

      treatmatchperavg[r, c] = Vector{Union{Float64, Missing}}(undef, lK);
      matchmatchperavg[r, c] = Vector{Union{Float64, Missing}}(undef, lK);
      treatoutcome[r, c] = Vector{Union{Float64, Missing}}(undef, lK);
      matchoutcome[r, c] = Vector{Union{Float64, Missing}}(undef, lK);
      possible[r, c] = [false for l in 1:lK];
      mdistes[r, c] = Vector{Union{Float64, Missing}}(undef, lK);
      [cals[cv][r, c] = Vector{Union{Float64, Missing}}(undef, lK) for cv in 1:length(cc.covariates)]

      for k in eachindex(1:lK)

        f = (k - 1) + cc.fmin;
        
        # need k, f, situation
        
        treatmatchperavg[r, c][k] = mean(tusdf[1:(trtdx + f - 1), cc.outcome]);
        matchmatchperavg[r, c][k] = mean(judat[1:(trtdx + f - 1), cc.outcome]);

        treatoutcome[r, c][k] = tusdf[trtdx + f, cc.outcome]
        matchoutcome[r, c][k] = judat[trtdx + f, cc.outcome]

        Lf = LF[k : ((k - 1) + lM)];
        c1j = judat[!, cc.t] .∈ Ref(Lf);
        # c2a = sum(c1j) >= min_match_len # minimum match period length
        c2 = !any(e -> e == 1, @view judat[c1j, cc.treatment]);

        if c2
          possible[r, c][k] = true;
    
          # do distance calculations
          if distances
    
            c1ousdf = @view judat[c1j, cc.covariates];
            c1i = tusdf[!, t] .∈ Ref(Lf);
            c1tusdf = @view tusdf[c1i, cc.covariates];
            
            mdist = [0.0]
            caldists = zeros(Float64, length(cc.covariates));
            lcnt = [0]
            
            distance_sums!(
              tmin, Σinvdict,
              cc.covariates, c1tusdf, c1ousdf, LF, mmin:mmax, caldists, mdist, lcnt
            )

            mdistes[r, c][k] = mdist[1] / lcnt[1];
            
            # add caliper distances

            for cnum = eachindex(covariates)
              # covate = covariates[cnum]
              cals[cnum][r,c][k] = (caldists[cnum] / lcnt[1]);
            end
          end
        end
      end
    end
  end

end




  # vect
  # @view vect[[1,2,3,3,2,1,]]

  L = length(cc.fmin:cc.fmax);

  vars = [
    :treatpretreat, :matchpretreat,
    :treatmatchperavg, :matchmatchperavg,
    :treatoutcome, :matchoutcome
  ];

  cc.matches[!, :possible] = [[false for i in 1:L] for i in 1:nrow(cc.matches)];

  for var in vars
    cc.matches[!, var] = [Vector{Union{Float64, Missing}}(undef, L) for i in 1:nrow(cc.matches)]
  end

  for tv in cc.covariates
    cc.matches[!, tv] = [Vector{Union{Float64, Missing}}(undef, L) for i in 1:nrow(cc.matches)]
  end

  cc.matches[!, :mdist] = [Inf .* Vector{Float64}(undef, L) for i in 1:nrow(cc.matches)]

  treatment_points = findall(dat[!, cc.treatment] .== 1);
  rid = unique(dat[treatment_points, cc.id]);

  sort!(cc.matches, [cc.t, cc.id, :matchunit])

  Σinvdict, σd = get_Σs(
      dat, rid, cc.covariates, cc.t, cc.id, cc.treatment
    );

  return cc, Σinvdict
end

function domatch!(cc::cicmodel, dat::DataFrame, Σinvdict; distances = true)
  lmm = length(cc.mmin:cc.mmax)

  gdf = groupby(cc.matches, [cc.t, cc.id]);
  tmin = minimum(dat[!, cc.t]);
  tref = cc.reference;

  @inbounds Threads.@threads for gi in gdf
  # for i = 1:10
  #  gi = @views gdf[i];
  
    iu = gi[1, cc.id];
    tt = gi[1, cc.t];

    ttl = tt + cc.fmin + cc.mmin;
    ttu = tt + cc.fmax;
    sdf = @view dat[(dat[!, cc.t] .>= ttl) .& (dat[!, cc.t] .<= ttu), :];
    tusdf = @view sdf[sdf[!, id] .== iu, :];

    # for @eachrow! cc.matches[parentindices(gi)[1], :] begin
    for r in 1:nrow(gi)
      judat = @view sdf[sdf[!, cc.id] .== gi.matchunit[r], :];
      # judat = sdf[sdf[!, cc.id] .== gi.matchunit[80], :];

      LF = (tt + cc.fmin) + cc.mmin : (tt + cc.fmax) + cc.mmax;
    
      tstrt = max((tt + cc.fmin) + cc.mmin, tmin) # begin of data or earliest L
      # this is the first time point in sdf
      # assumes same start date for all units
      trtdx = tt - tstrt + 1; # first index is 1

      if trtdx < 2 # treatment must be present, and day before treatment
        # possibles[ijk, :poss] = false; # default already
        continue
      end
    
      for k in eachindex(cc.fmin:cc.fmax)

        f = (k - 1) + cc.fmin;
        
        # need k, f, situation
        gi[r, :treatpretreat][k] = tusdf[trtdx + tref, cc.outcome];
        gi[r, :matchpretreat][k] = judat[trtdx + tref, cc.outcome];
  
        gi[r, :treatmatchperavg][k] = mean(tusdf[1:(trtdx + f - 1), cc.outcome]);
        gi[r, :matchmatchperavg][k] = mean(judat[1:(trtdx + f - 1), cc.outcome]);

        gi[r, :treatoutcome][k] = tusdf[trtdx + f, cc.outcome]
        gi[r, :matchoutcome][k] = judat[trtdx + f, cc.outcome]

        Lf = LF[k : ((k - 1) + lmm)];
        c1j = judat[!, cc.t] .∈ Ref(Lf);
        # c2a = sum(c1j) >= min_match_len # minimum match period length
        c2 = !any(e -> e == 1, judat[c1j, cc.treatment]);

        if c2
          gi[r, :possible][k] = true;
    
          # do distance calculations
          if distances
    
            c1ousdf = @view judat[c1j, cc.covariates];
            c1i = tusdf[!, t] .∈ Ref(Lf);
            c1tusdf = @view tusdf[c1i, cc.covariates];
            
            mdist = [0.0]
            caldists = zeros(Float64, length(cc.covariates));
            lcnt = [0]
            
            distance_sums!(
              tmin, Σinvdict,
              cc.covariates, c1tusdf, c1ousdf, LF, mmin:mmax, caldists, mdist, lcnt
            )

            gi[r, :mdist][k] = mdist[1] / lcnt[1];
            
            # add caliper distances

            for cnum = eachindex(covariates)
              covate = covariates[cnum]
              # covate = string(covate)
              gi[r, covate][k] = (caldists[cnum] / lcnt[1]);
            end
          end
        end
      end
    end

  end
  return cc
end

# 2858.395840 seconds (114.88 G allocations: 4.771 TiB, 50.77% gc time, 0.00% compilation time)
# @time domatch!(cc, dat, Σinvdict)

function matching!(cc, dat)
  @time cc.matches, cc, Σinvdict = setupmatching(cc, dat);

  # remove known to be unacceptable matches
  @time begin
    cc.matches[!, :include] .= false;
    for i in 1:nrow(cc.matches)
      if ceil(cc.matches[i, :fips], digits = -3) != ceil(cc.matches[i, :matchunit], digits = -3)
        cc.matches[i, :include] = true
      end
    end
    @time cc.matches = cc.matches[cc.matches.include, :];
    # cc.matches = @subset(cc.matches, :include .== true);
    # @select(cc.matches, Not(:include))
    @time domatch!(cc, dat, Σinvdict)
  end
return cc
end


@time matching!(cc, dat)

save_object("new_new_matches.jld2", cc)

keeps = Vector{Bool}(undef, nrow(cc.matches));
begin
    for i = 1:nrow(cc.matches)
        keeps[i] = any(cc.matches.possible[i])
    end
end

@view cc.matches[keeps, :mdist];

cc.covariates .∈ Ref(Symbol.(names(cc.matches)))


gg = gi[3054:end,:]

@eachrow gg begin
  for cnum = 1:3
    covate = (cc.covariates)[cnum]
    $(covate) .= 2;
  end
end
