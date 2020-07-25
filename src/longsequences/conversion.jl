###
### Conversion & Promotion
###
###
### Conversion methods for LongSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

###
### Promotion
###
for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{LongSequence{A}}, ::Type{LongSequence{B}}) where {A<:$alph,B<:$alph}
        return LongSequence{promote_rule(A, B)}
    end
end

###
### Conversion
###

function Base.convert(::Type{T}, seq::LongSequence{<:NucleicAcidAlphabet}) where
         {T<:LongSequence{<:NucleicAcidAlphabet}}
    return T(seq)
end

# Convert from a LongSequence to to a DNA or RNA vector
function Base.convert(::Type{LongSequence{A}}, seq::Vector) where A<:Alphabet
    return LongSequence{A}(seq)
end

Base.parse(::Type{LongSequence{A}}, seq::AbstractString) where A = LongSequence{A}(seq)
