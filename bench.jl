# benchmark the motif finding algorithms
include("utils.jl")
include("search.jl")
include("data.jl")


"""
Print benchmarks for motifsearch_greedy
"""
function bench_greedy()

    count = 16
    len = 300
    motiflen = 25

    scores = []
    times = []

    for i = 1:100
        seqs, motif, inds = gendata(count, len, motiflen, 2)
        truescore = score(seqs, inds, motiflen)
        ret = @timed score(seqs, motifsearch_greedy(seqs, motiflen), motiflen) / truescore
        push!(scores, ret[1])
        push!(times, ret[2])
        if i % 10 == 0
            @printf("Step %d\tMean time per execution: %gs\n", i, mean(times))
        end
    end

    @printf("Mean score for motifsearch_greedy: %g\n", mean(scores))

end
