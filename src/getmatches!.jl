# getmatches!.jl

function getmatches!(
  observations, matches,
  rg, trtg, ids, fmin, fmax; sliding = false
)
  
  for (ob, tob) in zip(observations, matches)
    (tt, tu) = ob;
    @unpack mus = tob;
    obsmatches!(mus, rg, trtg, ids, tt, tu, fmin, fmax; sliding = sliding);
  end
  return matches
end

function obsmatches!(mus, rg, trtg, uid, tt, tu, fmin, fmax; sliding = false)

  for (mu, rmus) in zip(uid, eachrow(mus))

    # testing
    # (mu, rmus) = collect(zip(uid, eachrow(mus)))[2000]

    if tu == mu
      rmus .= false
    else

      # require that panels are balanced
      # will need to handle missingness case eventually
      # data are padded
      
      # ftrue = @views tobs.mus[m, :];
      
      pollution = trtg[(tt, mu)];
      
      # testing for tt = 2
      # ftrue = fill(true, flen);
      # pollution = zeros(Int, 31);
      # pollution[findfirst(rg[(tt, mu)] .== 3)] = 1;

      if sum(pollution) > 0
        getfset!(
          rmus, fmin, fmax, pollution, rg[(tt, mu)], tt; sliding = sliding
        )

      # else -- leave as true
      #   tobs.mus[m] = true # leave
      #   tobs.fs[m] = ftrue

      end
    end
  end
  return mus
end

"""
    getfset!(ftrue, fmin, fmax, pollution, gt, tt; sliding = false)

The PO for a unit at t is 40 days before up to 10 days before.
This is the period in which an intervention can directly alter the death rate at t. The 41st day after an intervention is the first day in which we could observe *no* exogenous deaths. Where *no* means that we are above the 95th percentile of the infection-death distribution (as in days after intervention).

We are concerned with treatments whose outcome windows overlap the outcome window of a focal treatment event. For some day q+f after a treatment, we want to ensure that there are no treatment events in a match that take place 10 to 40 days before q+f.

e.g.

for f = fmin
  - we don't care about a treatment that occurs 31 days before a treatment
  - 40 - 10 + 1
for f = fmax
  - we don't care about a treatment that occurs 01 days before a treatment
  - 40 - 40 + 1 = fmax - f + 1

## old, worse description
All fs are allowable if the matchunit is treated no closer than 31 days (before or after) the treatment.
  => mu treatment at 31th day after treatment means first outcome
      window day is after tu's fmax
  => 41st day after mu treatment is beyond its outcome window
      so it is fully elgible to be a match to a unit with an fmin
      on that day. Then the treatment can happen 10 (fmin), days before.
"""
function getfset!(ftrue, fmin, fmax, pollution, gt, tt; sliding = false);

  # (d, τ) = collect(zip(pollution, gt))[10]
  for (d, τ) in zip(pollution, gt)
    # tt+fmin:tt+fmax, τ+fmin:τ+fmax
    if d > 0
      if τ == tt # cannot do
        for φ in eachindex(ftrue); ftrue[φ] = false; end
        return ftrue
      else
        if (tt - τ > 0) & (tt - τ < (fmax - fmin + 1))
          # if before focal treatment, and within 30 days

          # φ possibly over 1:31
          # not usable towards the end,
          # tt - τ gives the number of days before
          # that the treatment happens in the match
          # this is also the number of usable
          # outcome days at the end for the pairing
          # 1:(fmax - fmin + 1) - (tt - τ)

          # if sliding, reasoning above applies
          if sliding
            for φ in 1:(fmax - fmin + 1) - (tt - τ)
              ftrue[φ] = false
            end
          else
            # if not sliding, all are bad
            # since we work off fmin PO window
            for φ in eachindex(ftrue); ftrue[φ] = false end
          end


        elseif (tt - τ < 0) & (tt - τ < (fmin + fmax + 1))
          # if after focal treatment and within 30 days after
          # e.g. if there is a treatment 30 days after a treatment,
          # day 40 cannot be estimated with that match

          # φ possibly over 1:31
          # for post-treated-treatment
          # treatment in match unit,
          # these are not valid:
          # fmin + (τ - tt) - fmin + 1:(fmax - fmin + 1)

          # if sliding or not, since only later fs are affected
          # by later treatments in a potential match
          for φ in (τ - tt) + 1:(fmax - fmin + 1)
            ftrue[φ] = false
          end

        end
      end
    end
  end
  return ftrue
end
