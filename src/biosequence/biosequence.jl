abstract type BioSequence{A <: Alphabet} end

# Aliases and shorthands for describing subsets of the BioSequence type...
const NucleotideSeq = BioSequence{<:NucleicAcidAlphabet}


# Required traits and methods
# ---------------------------

# Base.length must be defined for every T<:BioSequence.

# As must the following...

"""
Return the data member of `seq` that stores the encoded sequence data.
"""
@inline function encoded_data(seq::BioSequence)
    error(
        string(
            "encoded_data has not been defined for BioSequence type: ",
            typeof(seq),
            ". It is required for any BioSequence subtype."
        )
    )
end


# Provided traits and methods
# ---------------------------

# These traits and methods are defined automatically for any subtype of BioSequence{A}.
# They may be overloaded for your concrete BioSequence sub-type if it is nessecery.

@inline encoded_data_eltype(seq::BioSequence) = eltype(encoded_data(seq))

"""
Return the `Alpahbet` type defining the possible biological symbols
and their encoding for a given biological sequence.
"""
@inline function Alphabet(::Type{<:BioSequence{A}}) where {A <: Alphabet}
    return A()
end

"Return the `Alpahbet` type that defines the biological symbols allowed for `seq`."
@inline function Alphabet(seq::BioSequence)
    return Alphabet(typeof(seq))
end

BioSymbols.alphabet(::Type{BioSequence{A}}) where {A<:Alphabet} = alphabet(A)

BitsPerSymbol(seq::BioSequence) = BitsPerSymbol(Alphabet(seq))
bits_per_symbol(seq::BioSequence) = bits_per_symbol(Alphabet(seq))

# The generic functions for any BioSequence...
include("indexing.jl")
include("conversion.jl")
include("predicates.jl")
include("operations.jl")
include("printing.jl")
include("transformations.jl")