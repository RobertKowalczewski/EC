include("helpers.jl")
include("moves.jl")
include("steepest-ls-candidates.jl")

using Random, Statistics

data_path = "lab1/TSPB.csv"
distance_matrix, costs = prepare_data(data_path)

# Benchmark: random starts; baseline vs candidates (same starts)
function benchmark_ls(distance_matrix, costs; runs=50, n_candidates=5, intra_fn=intra_two_edges_exchange)
    inter_fn = inter_two_nodes_exchange

    base_objs = Float64[]; base_times = Float64[]
    cand_objs = Float64[]; cand_times = Float64[]

    Random.seed!(42)
    for run in 1:runs
        start_sol = random_start(distance_matrix, costs, nothing, run)

        t = @elapsed begin
            obj_b, _ = steepest_local_search(start_sol, distance_matrix, costs, intra_fn, inter_fn)
            push!(base_objs, obj_b)
        end
        push!(base_times, t)

        t = @elapsed begin
            obj_c, _ = steepest_local_search_candidates(start_sol, distance_matrix, costs, intra_fn, inter_fn; n_candidates=n_candidates)
            push!(cand_objs, obj_c)
        end
        push!(cand_times, t)
    end

    fmt(x) = round(x; digits=5)

    println("Baseline (no candidates): Obj avg=$(fmt(mean(base_objs))) min=$(fmt(minimum(base_objs))) max=$(fmt(maximum(base_objs))) | Time[s] avg=$(fmt(mean(base_times)))")
    println("Candidates (k=$(n_candidates)): Obj avg=$(fmt(mean(cand_objs))) min=$(fmt(minimum(cand_objs))) max=$(fmt(maximum(cand_objs))) | Time[s] avg=$(fmt(mean(cand_times)))")
end

# Run benchmark; adjust n_candidates (default 10) experimentally
benchmark_ls(distance_matrix, costs; runs=200, n_candidates=3, intra_fn=intra_two_edges_exchange)
benchmark_ls(distance_matrix, costs; runs=200, n_candidates=5, intra_fn=intra_two_edges_exchange)
benchmark_ls(distance_matrix, costs; runs=200, n_candidates=10, intra_fn=intra_two_edges_exchange)
