# Utility functions for motif finding algorithms

DNA = Dict("A" => 1, "T" => 2, "C" => 3, "G" => 4)

function align(sequences, s, l)
    """
    Return the length(*s*) x *l* alignment matrix for *sequences*
    """
    align = []
    for i = 1:length(s)
        # add the motif substring to align array
        push!(align, sequences[i][s[i]:s[i]+l-1])
    end
    return align
end

function getprofile(alignment)
    """
    Return the profile matrix of an alignment
    """
    profile = zeros(Int32, 4, length(alignment[1]))
    for i = 1:length(alignment)
        for j = 1:length(alignment[1])
            # get the character at position [i,j] in alignment
            nucl = alignment[i][[j]]
            profile[DNA[nucl], j] += 1
        end
    end
    return profile
end

score(profile) = sum(maximum(profile, 1))

function score(sequences, s, l)
    """
    Find the alignment and calculate the score of the corresponding profile
    matrix.
    """
    alignment = align(sequences, s, l)
    profile = getprofile(alignment)
    return score(profile)
end
