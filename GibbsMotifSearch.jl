using Distributions
include("utils.jl")

function motifsearch_gibbs(sequences, l, max_steps_since)
    """
    Given a 2D array of strings *sequences*, randomized greedy 
    search for a motif of length *l*. *max_steps_since* indicates how many steps
    are taken after score stops increasing.
    Return an array of indices corresponding to the alignment found with the 
    highest score.
    """
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
    while steps_since_inc < max_steps_since
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
