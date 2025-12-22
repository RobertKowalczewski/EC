# Evolution from lab6 to lab10: ILS with Memory

## Overview of Changes

### **lab6: Basic ILS**
- Random initial solutions
- Full neighborhood evaluation in local search
- Simple perturbation (2-opt/swap/replace mix)
- Strict hill-climbing acceptance
- Each LS iteration evaluates ALL moves

### **lab10: ILS with Memory**
- Greedy initial solutions
- Candidate list restriction (k-nearest neighbors)
- Double-bridge perturbation
- Non-monotone acceptance (1% worse solution acceptance)
- **Memory-based incremental move evaluation**

---

## High-Level Algorithm Comparison

### **lab6 ILS (Basic)**
```
FUNCTION IteratedLocalSearch_lab6(distance_matrix, costs, time_limit):
    x ← RandomSolution()                           // Random start
    (objective, x) ← LocalSearch_Full(x)           // Full neighborhood
    
    WHILE time_remaining DO
        y ← SimplePerturb(x)                       // 2-opt/swap/replace
        (objective_y, y) ← LocalSearch_Full(y)     // Recalculate ALL moves
        
        IF objective_y < objective THEN            // Strict acceptance
            objective, x ← objective_y, y
        END IF
    END WHILE
    
    RETURN (objective, x)
END FUNCTION
```

### **lab10 ILS with Memory**
```
FUNCTION IteratedLocalSearch_lab10(distance_matrix, costs, time_limit, neighbors):
    x ← GreedySolution()                           // ✓ Better start
    best_x ← x
    best_obj ← ∞
    
    (objective, x) ← LocalSearch_Memory(x, neighbors)  // ✓ Candidate + Memory
    best_x ← x
    best_obj ← objective
    
    WHILE time_remaining DO
        y ← DoubleBridgePerturb(x)                // ✓ Stronger perturbation
        IF random() < 0.5 THEN
            y ← SmallNoise(y)                      // ✓ Additional diversification
        END IF
        
        (objective_y, y) ← LocalSearch_Memory(y, neighbors)  // ✓ Incremental updates
        
        IF objective_y < objective THEN
            objective, x ← objective_y, y
            IF objective < best_obj THEN          // ✓ Track best-ever
                best_obj, best_x ← objective, x
            END IF
        ELSE IF random() < 0.01 THEN              // ✓ Non-monotone acceptance
            objective, x ← objective_y, y
        END IF
    END WHILE
    
    RETURN (best_obj, best_x)
END FUNCTION
```

---

## Key Innovation: Memory-Based Local Search

### **lab6 Approach: Recalculate Everything**
```
FUNCTION LocalSearch_lab6(solution):
    improved ← true
    WHILE improved DO
        best_delta ← 0
        best_move ← null
        
        // Evaluate ALL O(n²) intra-route moves
        FOR i in 1..n, j in i+2..n DO
            delta ← evaluate_intra(i, j)          // Full calculation
            IF delta < best_delta THEN
                best_delta ← delta
                best_move ← (i, j)
            END IF
        END FOR
        
        // Evaluate ALL O(n·N) inter-route moves  
        FOR i in 1..n, node in available_nodes DO
            delta ← evaluate_inter(i, node)       // Full calculation
            IF delta < best_delta THEN
                best_delta ← delta
                best_move ← (i, node)
            END IF
        END FOR
        
        IF best_move exists THEN
            apply(best_move)
            improved ← true
        ELSE
            improved ← false
        END IF
    END WHILE
END FUNCTION

// Each iteration: O(n²) + O(n·N) evaluations
// Total: many iterations × full neighborhood evaluation
```

### **lab10 Approach: Cache + Incremental Update**
```
FUNCTION LocalSearch_lab10_Memory(solution, neighbors):
    
    // ✓ PHASE 1: PREPROCESSING (once per local search call)
    Precompute k=15 nearest neighbors for each node
    
    // ✓ PHASE 2: INITIALIZATION
    memory ← empty cache
    FOR each node u in solution DO
        FOR each neighbor v in neighbors[u] (k neighbors only) DO
            delta ← evaluate_move(u, v)
            memory[u,v] ← delta                   // Store in cache
        END FOR
    END FOR
    
    improved ← true
    WHILE improved DO
        // ✓ PHASE 3: LOOKUP (very fast)
        best_delta ← 0
        best_move ← null
        
        FOR each cached entry in memory DO       // Just O(n·k) lookups
            IF memory[entry] < best_delta THEN
                best_delta ← memory[entry]
                best_move ← entry
            END IF
        END FOR
        
        IF best_move exists THEN
            apply(best_move)
            
            // ✓ PHASE 4: INCREMENTAL UPDATE (only affected moves)
            affected_nodes ← identify_affected_by(best_move)  // Small subset
            
            FOR each node u in affected_nodes DO
                FOR each neighbor v in neighbors[u] DO
                    delta ← evaluate_move(u, v)
                    memory[u,v] ← delta          // Update only affected
                END FOR
            END FOR
            
            improved ← true
        ELSE
            improved ← false
        END IF
    END WHILE
END FUNCTION

// Initialization: O(n·k) evaluations where k=15
// Each iteration: 
//   - Lookup: O(n·k) comparisons (cached values)
//   - Update: O(affected·k) re-evaluations (affected << n)
// Speedup: ~10-50x
```

---

## Detailed Change Breakdown

### **1. Initial Solution Construction**

**lab6:**
```
x ← RandomPermutation()[1..n/2]
```

**lab10:**
```
FUNCTION GreedySolution():
    solution ← [random_start_node]
    WHILE size < n/2 DO
        Find node and position with minimum insertion cost
        Insert node at best position
    END WHILE
    RETURN solution
END FUNCTION
```
**Impact:** Much better starting solutions → fewer LS iterations needed

---

### **2. Neighborhood Restriction**

**lab6:**
```
// Evaluate moves between ALL node pairs
FOR i, j in all_combinations(nodes) DO
    evaluate_move(i, j)
END FOR
```

**lab10:**
```
// Preprocess once
neighbors ← PrecomputeKNearest(distance_matrix, k=15)

// Only evaluate moves to nearest neighbors
FOR i in solution DO
    FOR j in neighbors[i] DO  // Only k=15 neighbors
        evaluate_move(i, j)
    END FOR
END FOR
```
**Impact:** Reduces neighborhood size by ~10-100x with minimal quality loss

---

### **3. Move Evaluation Strategy**

**lab6:**
```
// Every iteration, recalculate from scratch
FOR each possible move DO
    delta ← ComputeDelta(move)  // Fresh calculation
    IF delta < best THEN
        best ← delta
    END IF
END FOR
```

**lab10:**
```
// Initialize once
FOR each possible move DO
    memory[move] ← ComputeDelta(move)
END FOR

// Each iteration: lookup cached values
WHILE improved DO
    best_move ← FindMin(memory)  // Fast lookup
    
    // Only update affected entries
    affected_moves ← GetAffectedBy(best_move)
    FOR move in affected_moves DO
        memory[move] ← ComputeDelta(move)  // Selective update
    END FOR
END WHILE
```
**Impact:** Eliminates redundant calculations, ~10-50x speedup per LS call

---

### **4. Perturbation Mechanism**

**lab6:**
```
FUNCTION Perturb(solution):
    FOR _ in 1..few_steps DO
        random_choice among:
            - 2-opt (reverse segment)
            - Swap two nodes
            - Replace one node
    END FOR
END FUNCTION
```

**lab10:**
```
FUNCTION DoubleBridgePerturb(solution):
    Pick 4 random cut points: p1 < p2 < p3 < p4
    Reorder: [1..p1] + [p3..p4] + [p2..p3] + [p1..p2] + [p4..n]
END FUNCTION

FUNCTION EnhancedPerturb(solution):
    y ← DoubleBridgePerturb(solution)
    IF random() < 0.5 THEN
        y ← SmallRandomNoise(y)
    END IF
    RETURN y
END FUNCTION
```
**Impact:** Stronger diversification, escapes local optima better

---

### **5. Acceptance Criterion**

**lab6:**
```
IF objective_new < objective_current THEN
    current ← new
END IF
```
**Result:** Gets stuck in local optima

**lab10:**
```
IF objective_new < objective_current THEN
    current ← new
    IF objective_new < objective_best THEN
        best ← new
    END IF
ELSE IF random() < 0.01 THEN  // 1% acceptance of worse
    current ← new
END IF

RETURN best  // Not current!
```
**Impact:** Explores wider search space, maintains best solution found

---

## Performance Comparison

### **Time Complexity per ILS Iteration**

| Component | lab6 | lab10 (no memory) | lab10 (with memory) |
|-----------|------|-------------------|---------------------|
| **Perturbation** | O(n) | O(n) | O(n) |
| **LS Initialization** | - | - | O(n·k) = O(15n) |
| **LS per iteration** | O(n² + n·N) | O(n·k) = O(15n) | O(n·k) lookup + O(affected·k) update |
| **LS total** | Many × O(n²) | Many × O(n·k) | Init + Few × O(small) |

### **Expected Speedup**
- **Candidate lists alone:** 10-20x faster
- **Memory + candidates:** 20-50x faster
- **Better initial solutions:** 2-5x fewer iterations needed
- **Overall:** ~50-100x faster convergence

### **Solution Quality**
- **lab6:** Good, but limited by time/iterations
- **lab10:** Better or equal due to:
  - Better starting points
  - More iterations possible (faster LS)
  - Better diversification (double-bridge)
  - Escape mechanism (non-monotone acceptance)

---

## Summary: What Makes Memory Version Special

1. **Caching Strategy**: Store all O(n·k) candidate move deltas in memory
2. **Fast Lookup**: Finding best move is O(n·k) comparison, not O(n·k) evaluation
3. **Incremental Updates**: Only recalculate affected moves (typically 5-20 nodes out of 100)
4. **Locality Exploitation**: Most moves stay valid/similar after small changes
5. **Amortization**: Initialization cost paid once, benefits multiply over iterations

---

## Core Memory Mechanism: The Innovation

### **The Problem with Traditional Local Search**
```
Traditional LS (lab6): Evaluate ALL moves every iteration

Iteration 1: Calculate delta for 10,000 moves → Pick best → Apply
Iteration 2: Calculate delta for 10,000 moves → Pick best → Apply  
Iteration 3: Calculate delta for 10,000 moves → Pick best → Apply
...

Problem: 99% of moves unchanged between iterations, but we recalculate everything!
```

### **The Memory Solution (lab10)**
```
Memory-based LS: Calculate once, update incrementally

Initialization: Calculate delta for 1,500 moves (only k=15 neighbors) → Store in cache

Iteration 1: Lookup 1,500 cached values → Pick best → Apply
            Update only ~100 affected moves in cache
            
Iteration 2: Lookup 1,500 cached values → Pick best → Apply
            Update only ~100 affected moves in cache
            
Iteration 3: Lookup 1,500 cached values → Pick best → Apply
            Update only ~100 affected moves in cache
...

Benefit: Fast lookups + small incremental updates instead of full recalculation
```

### **What Gets Cached**
```
FOR each node u in solution:
    FOR each of its k=15 nearest neighbors v:
        memory[u,v] = {
            delta: improvement value of this move
            type: INTRA or INTER move type
            args: move parameters (indices, nodes)
        }
```

### **Incremental Update Logic**
```
When we apply move affecting positions [i-1, i, i+1, i+2]:

1. Mark affected nodes:
   affected = {solution[i-1], solution[i], solution[i+1], solution[i+2]}

2. Update only entries involving affected nodes:
   FOR u in solution:
       IF u in affected OR any_neighbor(u) in affected:
           Recalculate memory[u, *]
       ELSE:
           Keep cached values (they're still valid!)
```

**Key insight:** After reversing segment [5..10], moves involving nodes 1,2,3 are unchanged!

---

## Visualization: Memory Impact

### **Without Memory (lab6 approach)**
```
Each LS iteration:
├─ Evaluate all n² intra moves     [expensive]
├─ Evaluate all n·N inter moves    [expensive]  
├─ Pick best                        [cheap]
└─ Apply move                       [cheap]

Cost per iteration: O(n²) evaluations
```

### **With Memory (lab10 approach)**
```
Initialization:
└─ Evaluate n·k candidate moves    [one-time cost]

Each LS iteration:
├─ Find best in cached values       [O(n·k) comparisons - FAST]
├─ Apply move                        [cheap]
└─ Update ~affected·k entries       [small subset]

Cost per iteration: O(n·k) lookups + O(small_subset) updates
```

---

## Complete lab10 Memory-Based Local Search (High Level)

```
FUNCTION LocalSearch_WithMemory(solution, distance_matrix, costs, neighbors):
    
    // SETUP: Position tracking
    pos[node] ← position of node in solution (0 if not in solution)
    
    // MEMORY INITIALIZATION: Cache all candidate moves
    memory ← empty cache of size (total_nodes × k_neighbors)
    
    FOR each node u in solution DO
        FOR k = 1 to 15 DO  // k nearest neighbors
            v ← neighbors[u][k]
            
            IF v in solution THEN
                // INTRA move: 2-opt between u and v
                delta ← calculate_2opt_delta(u, v)
                memory[u,k] ← {delta, INTRA, (pos[u], pos[v])}
            ELSE
                // INTER move: replace adjacent node with v
                delta ← calculate_best_replacement_delta(u, v)
                memory[u,k] ← {delta, INTER, (position, v)}
            END IF
        END FOR
    END FOR
    
    // MAIN LOOP: Use cached values
    improved ← true
    WHILE improved DO
        improved ← false
        
        // FAST LOOKUP: Find best cached move
        best_delta ← 0
        best_move ← null
        
        FOR each entry (u,k) in memory DO
            IF memory[u,k].delta < best_delta THEN
                best_delta ← memory[u,k].delta
                best_move ← (u, k)
            END IF
        END FOR
        
        IF best_move is null THEN
            BREAK  // Local optimum reached
        END IF
        
        // APPLY MOVE
        (u, k) ← best_move
        apply(memory[u,k].type, memory[u,k].args)
        
        // INCREMENTAL UPDATE: Only affected entries
        affected ← identify_nodes_with_changed_edges()
        
        FOR each node u in solution DO
            FOR k = 1 to 15 DO
                v ← neighbors[u][k]
                
                // Only update if u or v was affected
                IF u in affected OR v in affected THEN
                    // Recalculate this entry
                    IF v in solution THEN
                        delta ← calculate_2opt_delta(u, v)
                        memory[u,k] ← {delta, INTRA, (pos[u], pos[v])}
                    ELSE
                        delta ← calculate_best_replacement_delta(u, v)
                        memory[u,k] ← {delta, INTER, (position, v)}
                    END IF
                END IF
                // Else: cached value still valid, skip
            END FOR
        END FOR
        
        improved ← true
    END WHILE
    
    RETURN solution
END FUNCTION
```
