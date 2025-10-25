using Random

function steepest_local_search(initial_solution::TSPSolution, instance::TSPInstance, intra_type::Symbol; stats::SearchStats=SearchStats())
    current = copy_solution(initial_solution)
    while true
        stats.iterations += 1
        best_move = nothing
        best_delta = (0, 0)
        best_total = 0

        for move in generate_intra_moves(current, intra_type)
            delta = calculate_delta(current, instance, move)
            stats.evaluated_moves += 1
            total = delta[1] + delta[2]
            if total < best_total
                best_total = total
                best_move = move
                best_delta = delta
            end
        end

        for move in generate_inter_moves(current)
            delta = calculate_delta(current, instance, move)
            stats.evaluated_moves += 1
            total = delta[1] + delta[2]
            if total < best_total
                best_total = total
                best_move = move
                best_delta = delta
            end
        end

        if best_move === nothing
            break
        end

        apply_move!(current, instance, best_move, best_delta)
        stats.applied_moves += 1
    end

    return current, stats
end

function greedy_local_search(initial_solution::TSPSolution, instance::TSPInstance, intra_type::Symbol; rng::AbstractRNG=Random.default_rng(), stats::SearchStats=SearchStats())
    current = copy_solution(initial_solution)
    while true
        stats.iterations += 1
        categories = Symbol[:intra, :inter]
        shuffle!(rng, categories)
        improved = false

        for category in categories
            moves = category == :intra ? generate_intra_moves(current, intra_type) : generate_inter_moves(current)
            shuffle!(rng, moves)
            for move in moves
                delta = calculate_delta(current, instance, move)
                stats.evaluated_moves += 1
                total = delta[1] + delta[2]
                if total < 0
                    apply_move!(current, instance, move, delta)
                    stats.applied_moves += 1
                    improved = true
                    break
                end
            end
            if improved
                break
            end
        end

        if !improved
            break
        end
    end

    return current, stats
end