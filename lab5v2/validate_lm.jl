include("local-search.jl")
using Test

function make_saved_move(delta, orientations)
    SavedMove(delta, :intra,
        [(1, 2), (3, 4)],
        [(1, 3), (2, 4)],
        orientations)
end

solution1 = [1, 2, 3, 4]
pos1 = build_position_map(solution1)
removed_move = SavedMove(-1.0, :intra,
    [(1, 5)],
    [(1, 2)],
    [:forward])
@test check_if_applicable(solution1, pos1, removed_move)[1] === :remove

partial_solution = [2, 1, 3, 4]
pos_partial = build_position_map(partial_solution)
partial_move = make_saved_move(-2.0, [:forward, :forward])
@test check_if_applicable(partial_solution, pos_partial, partial_move)[1] === :skip

matching_solution = [1, 2, 3, 4]
pos_match = build_position_map(matching_solution)
match_move = make_saved_move(-3.0, [:forward, :forward])
status, new_solution = check_if_applicable(matching_solution, pos_match, match_move)
@test status === :apply
@test new_solution == [1, 3, 2, 4]

after_reversal = [2, 1, 4, 3]
pos_rev = build_position_map(after_reversal)
reverse_move = make_saved_move(-4.0, [:forward, :forward])
status_rev, new_solution_rev = check_if_applicable(after_reversal, pos_rev, reverse_move)
@test status_rev === :apply
@test new_solution_rev == [2, 4, 1, 3]

println("LM validation passed")
