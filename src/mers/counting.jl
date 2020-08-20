###
### Mer specific specializations of src/biosequence/counting.jl
###

for i in [(:_count_a, :count_a), (:_count_c, :count_c), (:_count_g, :count_g), (:_count_t, :count_t)]
    @eval begin
        @inline function $(i[1])(head::UInt64, tail...)
            return $(i[2])(head) + $(i[1])(tail...)
        end
        @inline $(i[1])() = 0
    end
end

@inline function _count_gc(head::UInt64, tail...)
    return gc_bitcount(head, BitsPerSymbol{2}()) + _count_gc(tail...)
end
@inline _count_gc() = 0

count_a(x::Kmer) = _count_a(x.data...) - n_unused(x)
count_c(x::Kmer) = _count_c(x.data...)
count_g(x::Kmer) = _count_g(x.data...)
count_t(x::Kmer) = _count_t(x.data...)
count_gc(x::Kmer) = _count_gc(x.data...)
Base.count(::typeof(isGC), x::Kmer) = count_gc(x)

@inline function Base.count(::typeof(!=), a::Kmer{A,K,N}, b::Kmer{A,K,N}) where {A,K,N}
    ad = a.data
    bd = b.data
    bits = ntuple(Val{N}()) do i
        Base.@_inline_meta
        @inbounds x = ad[i] ⊻ bd[i]
        return x
    end
    sum = 0
    @inbounds for i in 1:N
        sum += count_nonzero_bitpairs(bits[i])
    end
    return sum
end

@inline function Base.count(::typeof(==), a::Kmer{A,K,N}, b::Kmer{A,K,N}) where {A,K,N}
    ad = a.data
    bd = b.data
    bits = ntuple(Val{N}()) do i
        Base.@_inline_meta
        @inbounds x = ad[i] ⊻ bd[i]
        return x
    end
    sum = 0
    @inbounds for i in 1:N
        sum += count_00_bitpairs(bits[i])
    end
    return sum - n_unused(a)
end