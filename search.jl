using Distributions
include("utils.jl")

"""
Given a 2D array of strings *sequences*, greedy 
search for a motif of length *l*.
Return an array of indices corresponding to the 
alignment found with the highest score.
"""
function motifsearch_greedy(sequences, l)
    # TODO: Check that l <= sequence length
    bestmotif = ones(Int32, 2)
    
    for s₁ = 1:length(sequences[1])-l+1
        for s₂ = 1:length(sequences[2])-l+1
            # find best partial alignment of first 2 sequences 
            if score(sequences, [s₁, s₂], l) > score(sequences, bestmotif, l)
                bestmotif = [s₁, s₂]
            end
        end
    end
    s = copy(bestmotif)
    
    #iterate through the rest of the sequences
    for i = 3:length(sequences)
        # extend the default best motif
        push!(bestmotif, 1)
        for sᵢ = 1:length(sequences[i])-l+1
            if score(sequences, [s..., sᵢ], l) > score(sequences, bestmotif, l)
                # update the best motif
                bestmotif[i] = sᵢ
            end
        end
        push!(s, bestmotif[i])
    end
    return bestmotif
end


"""
Given a 2D array of strings *sequences*, randomized greedy 
search for a motif of length *l*.
Return an array of indices corresponding to the 
alignment found with the highest score.
"""
function motifsearch_random(sequences, l)
    # TODO: Check that l <= sequence length
    
    n = length(sequences[1])
    
    # generate an initial random alignment
    s = rand(1:n-l+1, length(sequences))
    
    # form profile P from s
    P = getprofile(align(sequences, s, l))
    
    bestscore = 0
    newscore = score(P)
    while newscore > bestscore
        bestscore = newscore
        for i = 1:length(sequences)
            # find a P-most probable l-mer from the ith sequence
            a = map(x -> probscore(x, P), [sequences[i][j:j+l-1] for j=1:n-l+1])
            # update the starting index for sᵢ
            s[i] = indmax(a)
        end
        # recalculate the profile
        P = getprofile(align(sequences, s, l))
        newscore = score(P)
    end
    
    return s
end


"""
Given a 2D array of strings *sequences*, randomized greedy 
search for a motif of length *l*. *leniency* indicates how many steps
are taken after score stops increasing.
Return an array of indices corresponding to the alignment found with the 
highest score.
"""
function motifsearch_gibbs(sequences, l, leniency)
    # TODO: Check that l <= sequence length
    
    n = length(sequences[1])
    
    # generate an initial random alignment
    s⃗ = rand(1:n-l+1, length(sequences))
    bests⃗ = copy(s⃗)

    bestscore = score(sequences, s⃗, l)
    currscore = bestscore
    # watch for convergence
    steps_since_inc = 0

    # do until convergence:
    while steps_since_inc < leniency
        # Randomly select one sequence, seqᵢ
        seqᵢ = rand(1:length(sequences))
        # Get the indices of the sequences left over
        inds = filter(x -> x != seqᵢ, 1:length(sequences))
        # Create profile from sequences not selected
        alignment = align(sequences[inds], s⃗[inds], l)
        profile = getprofile(alignment)
        # Calculate P-prob for each motif in seqᵢ
        a = map(x -> probscore(x, profile), [sequences[seqᵢ][j:j+l-1] for j=1:n-l+1])
        # if all zeros, set to 1 for uniform dist
        if all(a .== 0)
            a = ones(length(a))
        end
        p⃗ = a ./ sum(a)
        # Choose new starting position, sampled from probs found for seqᵢ
        c = Categorical(p⃗)
        s⃗[seqᵢ] = rand(c)

        # calculate new alignment score and check for convergence
        currscore = score(sequences, s⃗, l)
        if currscore <= bestscore
            # increment convergence check
            steps_since_inc += 1
        else 
            # update bests⃗ if higher score found
            bestscore = currscore
            bests⃗ = copy(s⃗)
            steps_since_inc = 0
        end

    end
    
    return bests⃗
end
