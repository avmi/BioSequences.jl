global reps = 10

@testset "Construction and Conversions" begin
    @test Kmer(DNA_A, DNA_G, DNA_T) === Kmer("AGT")
    @test Kmer(RNA_A, RNA_G, RNA_U) === Kmer("AGU")
    
    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(::Type{T}, seq::AbstractString) where {T<:Kmer}
        return String(T(seq)) == uppercase(seq)
    end

    # Check that DNAKmers can be constructed from a LongDNASeq
    #   LongDNASeq → Kmer → LongDNASeq
    function check_dnasequence_construction(::Type{T}, seq::LongDNASeq) where {T<:Kmer}
        return LongDNASeq(T(seq)) == seq
    end

    # Check that RNAKmers can be constructed from a LongRNASeq
    #   LongRNASeq → Kmer → LongRNASeq
    function check_rnasequence_construction(::Type{T}, seq::LongRNASeq) where {T<:Kmer}
        return LongRNASeq(T(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(::Type{T}, seq::LongSequence) where {T<:Kmer}
        return LongSequence(T(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
    function check_nucarray_kmer(::Type{M}, seq::Vector{T}) where {T<:NucleicAcid,M<:Kmer}
        return String([convert(Char, c) for c in seq]) == String(M(seq))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(::Type{T}, A2, seq::AbstractString) where {T<:Kmer}
        return String(LongSequence{A2}(T(LongSequence{A2}(seq)))) == uppercase(seq)
    end
    
    #=
    function check_uint_conversion(::Type{T}) where {T<:Kmer}
        U = BioSequences.encoded_data_type(T)
        uint = rand(typemin(U):U(one(U) << 2BioSequences.ksize(T) - 1))
        return convert(U, T(uint)) === uint
    end
    =#

    @testset "Kmer conversion" begin
        for len in [1, 16, 32, 64]
            
            if len <= 32
                # UInt64 conversions
                #@test all(Bool[check_uint_conversion(DNAKmer{len}) for _ in 1:reps])
                #@test all(Bool[check_uint_conversion(RNAKmer{len}) for _ in 1:reps])
                # String construction
                @test all(Bool[check_string_construction(DNAKmer{len}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_string_construction(RNAKmer{len}, random_rna_kmer(len)) for _ in 1:reps])
                # DNA/LongRNASeq Constructions
                @test all(Bool[check_dnasequence_construction(DNAKmer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
                @test all(Bool[check_rnasequence_construction(RNAKmer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
                # BioSequence Construction
                @test all(Bool[check_biosequence_construction(DNAKmer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
                @test all(Bool[check_biosequence_construction(RNAKmer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
                # Construction from nucleotide arrays
                @test all(Bool[check_nucarray_kmer(DNAKmer{len}, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
                @test all(Bool[check_nucarray_kmer(RNAKmer{len}, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
                # Roundabout conversions
                @test all(Bool[check_roundabout_construction(DNAKmer{len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(DNAKmer{len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(RNAKmer{len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(RNAKmer{len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            end
        end
    end

    @test_throws MethodError Kmer() # can't construct 0-mer using `Kmer()`
    @test_throws ArgumentError DNAKmer(dna"") # 0-mers not allowed
    @test_throws ArgumentError DNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws BioSequences.EncodeError Kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws BioSequences.EncodeError Kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws BioSequences.EncodeError RNAKmer(rna"ACGNU")# no Ns in kmers
    @test_throws BioSequences.EncodeError DNAKmer(dna"ACGNT") # no Ns in kmers
    @test_throws MethodError Kmer(RNA_A, DNA_A) # no mixing of RNA and DNA

    @testset "From strings" begin
        @test DNAKmer("ACTG") == DNAKmer(LongDNASeq("ACTG"))
        @test RNAKmer("ACUG") == RNAKmer(LongRNASeq("ACUG"))

        # N is not allowed in Kmers
        @test_throws Exception DNAMmer("ACGTNACGT")
        @test_throws Exception RNAKmer("ACGUNACGU")

        # Test string literals
        @test mer"ACTG"dna == DNAKmer(LongDNASeq("ACTG"))
        @test isa(mer"ACGT"dna, DNAKmer{4})
        @test_throws LoadError eval(:(mer"ACGN"dna))
        @test_throws LoadError eval(:(mer"ACG-"dna))
    end
    
    @testset "Capacity" begin
        @test BioSequences.capacity(DNAKmer(random_dna_kmer(10))) == 32
        @test BioSequences.capacity(RNAKmer(random_rna_kmer(10))) == 32
        @test BioSequences.capacity(DNAKmer(random_dna_kmer(32))) == 32
        @test BioSequences.capacity(RNAKmer(random_rna_kmer(32))) == 32
    end
    
    @testset "N unused" begin
        @test BioSequences.n_unused(DNAKmer(random_dna_kmer(10))) == 22
        @test BioSequences.n_unused(RNAKmer(random_rna_kmer(10))) == 22
        @test BioSequences.n_unused(DNAKmer(random_dna_kmer(32))) == 0
        @test BioSequences.n_unused(RNAKmer(random_rna_kmer(32))) == 0
    end
end
