# data generation functions for motif finding

DNA_arr = ["A", "T", "C", "G"]

# generate a random DNA sequence of length n
genseq(n) = join([DNA_arr[i] for i in rand(1:4, 10)])

function gendata(count, len, motiflen, num_mutations)
    """
    Generate *count* sequences of length *len*, each with a motif of length
    *motiflen* is inserted into a random index. *num_mutations* mutations are 
    added randomly to each motif in each sequence.
    
    Return an array of sequences, the motif, and the indices of the inserted 
    motifs in each sequence.
    """

    motif_arr = [DNA_arr[i] for i in rand(1:4, motiflen)]
    seqs_arr = [DNA_arr[i] for i in rand(1:4, count, len)]

    motif_ind = rand(1:len - motiflen + 1, count)

    # insert (possibly mutated) motifs into sequences
    for i = 1:count
        if num_mutations == 0
            seqs_arr[i, motif_ind[i]:motif_ind[i]+motiflen-1] = motif_arr
        else
            motif_arrᵢ = copy(motif_arr)
            for d = 1:num_mutations
                # add a random mutation to motif
                motif_arrᵢ[rand(1:motiflen)] = DNA_arr[rand(1:4)]
                
            end
            seqs_arr[i, motif_ind[i]:motif_ind[i]+motiflen-1] = motif_arrᵢ
        end
    end

    return [join(seqs_arr[i, :]) for i in 1:count], join(motif_arr), motif_ind
end
