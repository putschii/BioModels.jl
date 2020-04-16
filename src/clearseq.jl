function clearseq(Sequence)
        for i in 1:length(Sequence)
            pop!(Sequence)
        end
end
