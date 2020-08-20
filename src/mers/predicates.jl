###
### Mer specific specializations of src/biosequence/predicates.jl
###

Base.cmp(x::T, y::T) where {T<:Kmer} = cmp(x.data, y.data)
Base.:(==)(x::T, y::T) where {T<:Kmer} = x.data == y.data
Base.isless(x::T, y::T) where {T<:Kmer} = isless(x.data, y.data)

function Base.hash(x::Kmer{A,K,N}, h::UInt) where {A<:NucleicAcidAlphabet{2},K,N}
    #TODO: Plz review this as I honestly don't know about the behaviour of this hash:
    # does it have some specific desired properties regarding the values computed by the
    # hash function? My punt at a solution is just to do the operations over every word.
    # I believe this is the right thing to do as I reason this should basically produce
    # the same behaviour (assembly.... at least on my MBP) when N = 2 as when we used
    # to have UInt128 based big kmers.
    #return Base.hash(packed_data(x) ⊻ K ⊻ h)
    return Base.hash(map(xx -> xx ⊻ K ⊻ h, x.data))
end
