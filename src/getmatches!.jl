# getmatches!.jl

function getmatches!(
  observations, matches,
  rg, trtg, ids, fmin, fmax
)
  for (ob, tob) in zip(observations, matches)
    (tt, tu) = ob;
    @unpack mus = tob;
    obsmatches!(mus, rg, trtg, ids, tt, tu, fmin, fmax);
  end
  return matches
end

function obsmatches!(mus, rg, trtg, uid, tt, tu, fmin, fmax)

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
          rmus, fmin, fmax, pollution, rg[(tt, mu)], tt
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
      getfset!(ftrue, fmin, fmax, pollution, gt, tt)

All fs are allowable if the matchunit is treated no closer than 31 days (before or after) the treatment.
  => mu treatment at 31th day after treatment means first outcome
      window day is after tu's fmax
  => 41st day after mu treatment is beyond its outcome window
      so it is fully elgible to be a match to a unit with an fmin
      on that day. Then the treatment can happen 10 (fmin), days # #    before.
"""
function getfset!(ftrue, fmin, fmax, pollution, gt, tt);

  # (d, τ) = collect(zip(pollution, gt))[10]

  for (d, τ) in zip(pollution, gt)
    # tt+fmin:tt+fmax, τ+fmin:τ+fmax
    if d > 0
      if τ == tt
        for φ in eachindex(ftrue); ftrue[φ] = false; end
        return ftrue
      else
        if (tt - τ > 0) & (tt - τ < (fmax - fmin + 1))

          # φ possibly over 1:31
          # not usable towards the end,
          # tt - τ gives the number of days before
          # that the treatment happens in the match
          # this is also the number of usable
          # outcome days at the end for the pairing
          # 1:(fmax - fmin + 1) - (tt - τ)
          for φ in 1:(fmax - fmin + 1) - (tt - τ)
            ftrue[φ] = false
          end
        elseif (tt - τ < 0) & (tt - τ < 31)

          # φ possibly over 1:31
          # for post-treated-treatment
          # treatment in match unit,
          # these are not valid:
          # fmin + (τ - tt) - fmin + 1:(fmax - fmin + 1)
          for φ in (τ - tt) + 1:(fmax - fmin + 1)
            ftrue[φ] = false
          end
        end
      end
    end
  end
  return ftrue
end
