#!/usr/bin/env julia

"""
Run all validation steps in sequence and print a compact summary.

Steps:
1) Seed sweep (report-only by default; optional gates via --gate/TSCS_GATES)
2) Coverage (null DGP) with gate on overall coverage in [0.93, 0.97]
3) Placebo/permutation (real data) with gate on overall Type I ∈ [0.03, 0.07]

Usage examples:
  julia --project=. scripts/run_validations.jl
  julia --project=. scripts/run_validations.jl --out-dir test/validation --seeds 20 --iterations 400 --permutations 300
  TSCS_GATES=1 julia --project=. scripts/run_validations.jl --gate 1 --mae-thresh 0.03 --mxe-thresh 0.085

"""

"""
Quick usage

- Default run:
    - julia --project=. scripts/run_validations.jl
- With options:
    - julia --project=. scripts/run_validations.jl --out-dir test/validation
--seeds 20 --iterations 400 --permutations 300
- With seed-sweep gates:
    - TSCS_GATES=1 julia --project=. scripts/run_validations.jl --gate 1
--mae-thresh 0.03 --mxe-thresh 0.085

This should make local and CI runs straightforward while keeping to the
standard practice you chose.
"""

using Printf
using Dates

function parse_args()
    args = Dict{String,String}()
    i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            key = ARGS[i][3:end]
            i += 1
            i <= length(ARGS) || error("Missing value for --$key")
            args[key] = ARGS[i]
        else
            error("Unexpected arg: " * ARGS[i])
        end
        i += 1
    end
    return args
end

function ensure_dir(p::AbstractString)
    isdir(p) || mkpath(p)
    return p
end

function run_cmd(cmd::Cmd)
    t0 = time()
    try
        run(cmd)
        return (success=true, secs=time()-t0)
    catch e
        @error "Command failed" cmd string(e)
        return (success=false, secs=time()-t0)
    end
end

function extract_metric(path::AbstractString, key::AbstractString)
    # Minimal metric extraction from JSON-like output using regex; robust to spacing
    if !isfile(path)
        return nothing
    end
    s = read(path, String)
    m = match(Regex("\"$(key)\"\\s*:\\s*([-+0-9.eE]+)"), s)
    return m === nothing ? nothing : parse(Float64, m.captures[1])
end

function main()
    a = parse_args()
    out_dir = get(a, "out-dir", "test/validation") |> ensure_dir
    seeds = parse(Int, get(a, "seeds", "20"))
    iterations = parse(Int, get(a, "iterations", "400"))
    permutations = parse(Int, get(a, "permutations", "200"))

    # Optional seed-sweep gates
    gate = get(a, "gate", get(ENV, "TSCS_GATES", "0")) == "1"
    mae_thresh = get(a, "mae-thresh", "0.03")
    mxe_thresh = get(a, "mxe-thresh", "0.06")

    seed_out = joinpath(out_dir, "seed_sweep.json")
    cov_out = joinpath(out_dir, "bias_coverage.json")
    plc_out = joinpath(out_dir, "placebo.json")

    println("[run_validations] Output dir: $out_dir")

    # 1) Seed sweep (report-only by default)
    seed_cmd_parts = [
        Base.julia_cmd().exec[1], "--project=.", "scripts/seed_sweep.jl",
        "--seeds", string(seeds),
        "--out", seed_out,
    ]
    if gate
        append!(seed_cmd_parts, ["--gate", "1", "--mae-thresh", mae_thresh, "--mxe-thresh", mxe_thresh])
    end
    println("[run_validations] Running seed sweep…")
    r1 = run_cmd(Cmd(seed_cmd_parts))

    # 2) Coverage (null DGP) — gated
    cov_cmd = Cmd([
        Base.julia_cmd().exec[1], "--project=.", "scripts/sim_bias_coverage.jl",
        "--seeds", string(seeds), "--iterations", string(iterations),
        "--out", cov_out
    ])
    println("[run_validations] Running bias+coverage…")
    r2 = run_cmd(cov_cmd)

    # 3) Placebo (real data) — gated
    plc_cmd = Cmd([
        Base.julia_cmd().exec[1], "--project=.", "scripts/placebo_permutation.jl",
        "--permutations", string(permutations), "--iterations", string(iterations),
        "--out", plc_out
    ])
    println("[run_validations] Running placebo/permutation…")
    r3 = run_cmd(plc_cmd)

    # Summaries
    mae_mean = extract_metric(seed_out, "mae_mean")
    mxe_max = extract_metric(seed_out, "mxe_max")
    coverage = extract_metric(cov_out, "overall_coverage")
    type1 = extract_metric(plc_out, "overall_type1")

    println("\n================ Validation Summary ================")
    @printf("Seed sweep:   mae_mean=%s  mxe_max=%s  (%0.1fs)\n",
        isnothing(mae_mean) ? "NA" : @sprintf("%.4f", mae_mean),
        isnothing(mxe_max) ? "NA" : @sprintf("%.4f", mxe_max), r1.secs)
    @printf("Coverage:     overall=%s          (%0.1fs)\n",
        isnothing(coverage) ? "NA" : @sprintf("%.4f", coverage), r2.secs)
    @printf("Placebo T1:   overall=%s          (%0.1fs)\n",
        isnothing(type1) ? "NA" : @sprintf("%.4f", type1), r3.secs)
    println("Artifacts:")
    println("  • ", seed_out)
    println("  • ", cov_out)
    println("  • ", plc_out)
    println("===================================================\n")

    # Aggregate exit code: fail if any gated step failed
    exit_code = (r2.success && r3.success) ? 0 : 1
    if exit_code == 0
        println("[run_validations] All gated validations passed.")
    else
        println("[run_validations] One or more gated validations failed.")
    end
    exit(exit_code)
end

abspath(PROGRAM_FILE) == @__FILE__ && main()

