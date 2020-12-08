@inline encoded_data(x::Kmer) = x.data

@inline function bitindex(x::Kmer, i)
    i′ = i + div(64N - 2K, 2)
end