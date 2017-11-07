include("utils.jl")

function motifsearch_gibbs(sequences, l)
    """
    Given a 2D array of strings *sequences*, randomized greedy 
    search for a motif of length *l*.
    Return an array of indices corresponding to the 
    alignment found with the highest score.
    """
    # TODO: Check that l <= sequence length
    
    n = length(sequences[1])
    
    # generate an initial random alignment
    s = rand(1:n-l+1, length(sequences))

    # do until convergence:
        sᵢ = rand(1:length(sequences))
        # TODO: Randomly select one sequence, s' 
        # TODO: Create profile from sequences not selected
        # TODO: Calculate P-prob for each position in s'
        # TODO: Choose new starting position, sampled from probs found for s'
    
    

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
            # recalculate the profile
            P = getprofile(align(sequences, s, l))
            newscore = score(P)
        end
    end
    
    return s
end
