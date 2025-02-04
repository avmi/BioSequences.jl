@testset "LongSubSeq" begin

@testset "Construction" begin
	seq = LongSequence{AminoAcidAlphabet}([AA_A, AA_V, AA_W, AA_Y, AA_H])
	v1 = LongSubSeq{AminoAcidAlphabet}(seq.data, 2:4)
	v2 = LongSubSeq(seq, 2:4)
	v3 = view(seq, 2:4)
	v4 = @view seq[2:4]
	v5 = LongSubSeq(seq, :)
	vv = LongSubSeq(v1, 2:3)
	vv2 = v1[2:3]
	vv3 = LongSubSeq(vv)

	@test_throws BoundsError LongSubSeq(seq, 0:4)
	@test_throws BoundsError LongSubSeq(seq, 1:6)
	@test_throws BoundsError LongSubSeq(v1, 1:4)

    @test collect(v5) == collect(seq)
	@test typeof(v1) == typeof(v2) == typeof(v3) == typeof(v4) == typeof(vv) == LongSubSeq{AminoAcidAlphabet}
	@test v1 == v2 == v3 == v4

	@test vv == vv2 == vv3
	vv[1] = AA_V
	@test vv == vv2 == vv3

	@test LongSubSeq{AminoAcidAlphabet}(seq) == view(seq, eachindex(seq))

	s = dna"TASCTAWTA"
	@test collect(LongSubSeq{RNAAlphabet{4}}(s, 3:7)) == collect(LongRNA{4}(s)[3:7])
end

@testset "Basics" begin
	seq = LongSequence{AminoAcidAlphabet}([AA_A, AA_V, AA_W, AA_Y, AA_H])
	v1 = LongSubSeq{AminoAcidAlphabet}(seq, 2:4)
	v2 = LongSubSeq(seq, 1:0)

	@test length(v1) == 3
	@test !isempty(v1)
	@test isempty(v2)
	@test length(v2) == 0

	v1[1] = 'N'
	v1[2] = 'K'
	@test String(seq) == "ANKYH"
end

@testset "Equality" begin
    for size in [41, 504, 7]
        for offset in [0, 3, 32]
            seq = randrnaseq(size)
            seq2 = view(randrnaseq(offset) * seq * randrnaseq(15), offset+1:offset+size)
            @test seq == seq2
            seq[4] = RNA_Gap
            @test seq != seq2
            seq3 = view(seq, 1:size)
            @test seq == seq3
        end
    end
end

# Added after issue 260
@testset "Random construction" begin
	for i in 1:100
		seq = randdnaseq(rand(15:65))
		begin_ = min(lastindex(seq), rand(10:30))
		range = begin_:min(lastindex(seq), begin_ + rand(0:40))
		seq2 = view(seq, range)
		@test seq2 isa LongSubSeq{typeof(Alphabet(seq))}
		seq3 = LongSequence(seq2)
		@test typeof(seq) == typeof(seq3)
		@test seq[range] == seq2 == seq3
	end

	# See issue 260
	seq = dna"CATTTTTTTTTTTTTTT"
	seq2 = LongSequence(LongSubSeq(seq, 1:17))
	@test seq == seq2
end

@testset "Conversion" begin
	seq = LongDNA{4}("TAGTATCGAAMYCGNA")
	v = LongSubSeq(seq, 3:14)

	@test LongSequence(v) == seq[3:14]
	s2 = LongSequence{RNAAlphabet{4}}(seq[3:14])
	@test LongSequence{RNAAlphabet{4}}(v) == s2

	@test LongSequence(LongSubSeq{RNAAlphabet{4}}(seq)) == LongSequence{RNAAlphabet{4}}(seq)
	@test LongSubSeq{RNAAlphabet{4}}(seq) == LongSequence{RNAAlphabet{4}}(seq)
end

@testset "Transformations" begin
	# Reverse!
	str = "SKVANNSFDGRKIQAWPSRQ"
	seq = LongAA(str)
	seq2 = copy(seq)
	v = view(seq, 1:lastindex(seq))

	reverse!(v)
	@test seq == LongAA(reverse(str))
	@test seq == v
	@test seq != seq2

	reverse!(v)
	@test seq == LongAA(str)
	@test seq == v
	@test seq == seq2

	seq = LongDNA{4}("TGAGTCGTAGGAAGGACCTAAA")
	seq2 = copy(seq)
	v = LongSubSeq(seq2, 3:15)
	complement!(v)
	@test seq2[3:15] == complement(seq[3:15])
	@test seq2[1:2] == dna"TG"
	@test seq2[16:end] == seq[16:end]

	# A longer example to engage some inner loops
	seq = randdnaseq(38)
	seq2 = copy(seq)
	v = LongSubSeq(seq, 3:36)
	complement!(v)
	@test v == complement(seq2[3:36])
end

@testset "Copying" begin
	seq = LongRNA{4}("UAUUAACCGGAGAUCAUUCAGGUAA")
	v1 = view(seq, 1:3)
	v2 = view(seq, 4:11)

	@test copy(v1) == seq[1:3]
	@test copy(v2) == seq[4:11]

	# Can't resize views
	@test_throws Exception copy!(v1, v2)

	# Works even when sharing underlying data
	# Note it is not possible to have v2 == v3 after copying
	v3 = view(seq, 1:8)
	before = LongSequence(v3)
	copy!(v2, v3)
	@test v2 == before

	# Also works with nonshared data
	v1 = view(randaaseq(20), 3:17)
	v2 = view(randaaseq(20), 3:17)
	copy!(v1, v2)
	@test v1 == v2
end

end # seqview
