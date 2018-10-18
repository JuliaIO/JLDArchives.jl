using HDF5, JLD, Compat, Compat.Test, Compat.LinearAlgebra, LegacyStrings
using Compat: @warn

# Define variables of different types
x = 3.7
A = reshape(collect(1:15), 3, 5)
Aarray = Vector{Float64}[[1.2,1.3],[2.2,2.3,2.4]]
str = "Hello"
stringsA = String["It", "was", "a", "dark", "and", "stormy", "night"]
stringsU = String["It", "was", "a", "dark", "and", "stormy", "night"]
strings16 = convert(Array{UTF16String}, stringsA)
strings16_2d = reshape(strings16[1:6], (2,3))
empty_string = ""
empty_string_array = String[]
empty_array_of_strings = String[""]
tf = true
TF = A .> 10
B = [-1.5 sqrt(2) NaN 6;
     0.0  Inf eps() -Inf]
AB = Any[A, B]
t = (3, "cat")
c = ComplexF32(3,7)
cint = 1+im  # issue 108
C = Array(reinterpret(ComplexF64, vec(B)))
emptyA = zeros(0,2)
emptyB = zeros(2,0)
try
    global MyStruct
    mutable struct MyStruct
        len::Int
        data::Array{Float64}
        MyStruct(len::Int) = new(len)
        MyStruct(len::Int, data::Array{Float64}) = new(len, data)
    end
catch
end
ms = MyStruct(2, [3.2, -1.7])
msempty = MyStruct(5, Float64[])
sym = :TestSymbol
syms = [:a, :b]
d = Dict([(syms[1],"aardvark"), (syms[2], "banana")])
ex = quote
    function incrementby1(x::Int)
        x+1
    end
end
T = UInt8
char = 'x'
unicode_char = '\U10ffff'
α = 22
β = Any[[1, 2], [3, 4]]  # issue #93
vv = Vector{Int}[[1,2,3]]  # issue #123
typevar = Array{Int}[[1]]
typevar_lb = Vector{<:Integer}[[1]]
typevar_ub = (Vector{U} where U>:Int)[[1]]
typevar_lb_ub = (Vector{U} where Int<:U<:Real)[[1]]
# Immutable type:
rng = 1:5
# Type with a pointer field (#84)
struct ObjWithPointer
    a::Ptr{Cvoid}
end
objwithpointer = ObjWithPointer(0)
# Custom BitsType (#99)
primitive type MyBT 64 end
bt = reinterpret(MyBT, Int64(55))
# Symbol arrays (#100)
sa_asc = [:a, :b]
sa_utf8 = [:α, :β]
# SubArray (to test tuple type params)
subarray = view([1:5;], 1:5)
# Array of empty tuples (to test tuple type params)
arr_empty_tuple = Tuple{}[]
struct EmptyImmutable end
emptyimmutable = EmptyImmutable()
arr_emptyimmutable = [emptyimmutable]
mutable struct EmptyType end
emptytype = EmptyType()
arr_emptytype = [emptytype]
struct EmptyII
    x::EmptyImmutable
end
emptyii = EmptyII(EmptyImmutable())
struct EmptyIT
    x::EmptyType
end
emptyit = EmptyIT(EmptyType())
mutable struct EmptyTI
    x::EmptyImmutable
end
emptyti = EmptyTI(EmptyImmutable())
mutable struct EmptyTT
    x::EmptyType
end
emptytt = EmptyTT(EmptyType())
struct EmptyIIOtherField
    x::EmptyImmutable
    y::Float64
end
emptyiiotherfield = EmptyIIOtherField(EmptyImmutable(), 5.0)

# Unicode type field names (#118)
mutable struct MyUnicodeStruct☺{τ}
    α::τ
    ∂ₓα::τ
    MyUnicodeStruct☺{τ}(α::τ, ∂ₓα::τ) where {τ} = new{τ}(α, ∂ₓα)
end
unicodestruct☺ = MyUnicodeStruct☺{Float64}(1.0, -1.0)
# Arrays of matrices (#131)
array_of_matrices = Matrix{Int}[[1 2; 3 4], [5 6; 7 8]]
# Tuple of arrays and bitstype
tup = (1, 2, [1, 2], [1 2; 3 4], bt)
# Empty tuple
empty_tup = ()
# Non-pointer-free immutable
struct MyImmutable{T}
    x::Int
    y::Vector{T}
    z::Bool
end
nonpointerfree_immutable_1 = MyImmutable(1, [1., 2., 3.], false)
nonpointerfree_immutable_2 = MyImmutable(2, Any[3., 4., 5.], true)
struct MyImmutable2
    x::Vector{Int}
    MyImmutable2() = new()
end
nonpointerfree_immutable_3 = MyImmutable2()
# Immutable with a non-concrete datatype (issue #143)
struct Vague
    name::String
end
vague = Vague("foo")
# Immutable with a union of BitsTypes
struct BitsUnion
    x::Union{Int64, Float64}
end
bitsunion = BitsUnion(5.0)
# Immutable with a union of Types
struct TypeUnionField
    x::Union{Type{Int64}, Type{Float64}}
end
typeunionfield = TypeUnionField(Int64)
# Generic union type field
struct GenericUnionField
    x::Union{Vector{Int},Int}
end
genericunionfield = GenericUnionField(1)
# Array references
arr_contained = [1, 2, 3]
arr_ref = typeof(arr_contained)[]
push!(arr_ref, arr_contained, arr_contained)
# Object references
mutable struct ObjRefType
    x::ObjRefType
    y::ObjRefType
    ObjRefType() = new()
    ObjRefType(x, y) = new(x, y)
end
ref1 = ObjRefType()
obj_ref = ObjRefType(ObjRefType(ref1, ref1), ObjRefType(ref1, ref1))
# Immutable that requires padding between elements in array
struct PaddingTest
    x::Int64
    y::Int8
end
padding_test = PaddingTest[PaddingTest(i, i) for i = 1:8]
# Empty arrays of various types and sizes
empty_arr_1 = Int[]
empty_arr_2 = Array{Int}(undef, 56, 0)
empty_arr_3 = Any[]
empty_arr_4 = Array{Any}(undef, 0, 97)
# Moderately big dataset (which will be mmapped)
bigdata = [1:10000;]
# BigFloats and BigInts
bigints = big(3).^(1:100)
bigfloats = big(3.2).^(1:100)
# None
none = Union{}
nonearr = Array{Union{}}(undef, 5)
# nothing
scalar_nothing = nothing
vector_nothing = Union{Int,Nothing}[1,nothing]

# some data big enough to ensure that compression is used:
Abig = kron(Matrix(1.0I, 10, 10), reshape(1:400, 20, 20))
Bbig = Any[i for i=1:3000]
Sbig = "A test string "^1000

# Bitstype type parameters
mutable struct BitsParams{x}; end
bitsparamfloat  = BitsParams{1.0}()
bitsparambool   = BitsParams{true}()
bitsparamsymbol = BitsParams{:x}()
bitsparamint    = BitsParams{1}()
bitsparamuint   = BitsParams{0x01}()
bitsparamint16  = BitsParams{Int16(1)}()

# Tuple of tuples
tuple_of_tuples = (1, 2, (3, 4, [5, 6]), [7, 8])

iseq(x,y) = isequal(x,y)
iseq(x::MyStruct, y::MyStruct) = (x.len == y.len && x.data == y.data)
iseq(x::MyImmutable, y::MyImmutable) = (isequal(x.x, y.x) && isequal(x.y, y.y) && isequal(x.z, y.z))
@static if VERSION ≥ v"0.7.0-DEV.3693" # empty mutable structs are no longer singletons
    iseq(x::EmptyType, y::EmptyType) = true
    iseq(x::EmptyIT, y::EmptyIT) = true
    iseq(x::Array{EmptyType}, y::Array{EmptyType}) = size(x) == size(y)
    iseq(x::BitsParams{T}, y::BitsParams{T}) where {T} = true
    iseq(x::BitsParams, y::BitsParams) = false
end
iseq(x::Union{EmptyTI, EmptyTT}, y::Union{EmptyTI, EmptyTT}) = iseq(x.x, y.x)
iseq(c1::Array{Base.Sys.CPUinfo}, c2::Array{Base.Sys.CPUinfo}) = length(c1) == length(c2) && all([iseq(c1[i], c2[i]) for i = 1:length(c1)])
function iseq(c1::Base.Sys.CPUinfo, c2::Base.Sys.CPUinfo)
    for n in fieldnames(Base.Sys.CPUinfo)
        if getfield(c1, n) != getfield(c2, n)
            return false
        end
    end
    true
end
iseq(x::MyUnicodeStruct☺, y::MyUnicodeStruct☺) = (x.α == y.α && x.∂ₓα == y.∂ₓα)
iseq(x::Array{Union{}}, y::Array{Union{}}) = size(x) == size(y)
macro check(fid, sym)
    ex = quote
        let tmp
            try
                tmp = read($fid, $(string(sym)))
            catch e
                @warn("Error reading ", $(string(sym)))
                rethrow(e)
            end
            if !iseq(tmp, $sym)
                written = $sym
                error("For ", $(string(sym)), ", read value $tmp does not agree with written value $written")
            end
            written_type = typeof($sym)
            if typeof(tmp) != written_type
                error("For ", $(string(sym)), ", read type $(typeof(tmp)) does not agree with written type $(written_type)")
            end
        end
    end
    esc(ex)
end

# Test for equality of expressions, skipping line numbers
checkexpr(a, b) = @assert a == b
function checkexpr(a::Expr, b::Expr)
    @assert a.head == b.head
    i = 1
    j = 1
    while i <= length(a.args) && j <= length(b.args)
        if (isa(a.args[i], Expr) && a.args[i].head == :line) || isa(a.args[i], LineNumberNode)
            i += 1
            continue
        end
        if (isa(b.args[j], Expr) && b.args[j].head == :line) || isa(b.args[j], LineNumberNode)
            j += 1
            continue
        end
        checkexpr(a.args[i], b.args[j])
        i += 1
        j += 1
    end
    @assert i >= length(a.args) && j >= length(b.args)
end

for fname in ("v0.4.13-julia-0.3.jld", "v0.4.14-julia-0.4.0-dev+4483.jld")
    undefv = Array{Any}(undef, 1)
    undefm = Array{Any}(undef, 2, 2)
    ms_undef = MyStruct(0)

    fidr = jldopen(joinpath(dirname(@__FILE__), fname), "r")
    @check fidr x
    @check fidr A
    dsetA = fidr["A"]
    @test ndims(dsetA) == ndims(A)
    @test size(dsetA) == size(A)
    @test size(dsetA, 1) == size(A, 1)
    @test size(dsetA, 2) == size(A, 2)
    @test size(dsetA, 3) == size(A, 3)
    @check fidr Aarray
    @check fidr str
    @check fidr stringsA
    @check fidr stringsU
    @check fidr strings16
    @check fidr strings16_2d
    @check fidr empty_string
    @check fidr empty_string_array
    @check fidr empty_array_of_strings
    @check fidr tf
    @check fidr TF
    @check fidr AB
    @check fidr t
    @check fidr c
    @check fidr cint
    @check fidr C
    @check fidr emptyA
    @check fidr emptyB
    @check fidr ms
    @check fidr msempty
    @check fidr sym
    @check fidr syms
    @check fidr d
    exr = read(fidr, "ex")   # line numbers are stripped, don't expect equality
    checkexpr(ex, exr)
    @check fidr T
    @check fidr char
    @check fidr unicode_char
    @check fidr α
    @check fidr β
    @check fidr vv
    @check fidr rng
    @check fidr typevar
    @check fidr typevar_lb
    @check fidr typevar_ub
    @check fidr typevar_lb_ub

    # Special cases for reading undefs
    undefv = read(fidr, "undef")
    if !isa(undefv, Array{Any, 1}) || length(undefv) != 1 || isassigned(undefv, 1)
        error("For undef, read value does not agree with written value")
    end
    undefm = read(fidr, "undefs")
    if !isa(undefm, Array{Any, 2}) || length(undefm) != 4 || any(map(i->isassigned(undefm, i), 1:4))
        error("For undefs, read value does not agree with written value")
    end
    ms_undef = read(fidr, "ms_undef")
    if !isa(ms_undef, MyStruct) || ms_undef.len != 0 || isdefined(ms_undef, :data)
        error("For ms_undef, read value does not agree with written value")
    end

    @check fidr bt
    @check fidr sa_asc
    @check fidr sa_utf8
    # @check fidr subarray
    @check fidr arr_empty_tuple
    @check fidr emptyimmutable
    @check fidr emptytype
    @check fidr arr_emptyimmutable
    @check fidr arr_emptytype
    @check fidr emptyii
    @check fidr emptyit
    @check fidr emptyti
    @check fidr emptytt
    @check fidr emptyiiotherfield
    @check fidr unicodestruct☺
    @check fidr array_of_matrices
    @check fidr tup
    @check fidr empty_tup
    @check fidr nonpointerfree_immutable_1
    @check fidr nonpointerfree_immutable_2
    @check fidr nonpointerfree_immutable_3
    vaguer = read(fidr, "vague")
    @test typeof(vaguer) == typeof(vague) && vaguer.name == vague.name
    @check fidr bitsunion
    @check fidr typeunionfield
    @check fidr genericunionfield

    arr = read(fidr, "arr_ref")
    @test arr == arr_ref
    @test arr[1] === arr[2]

    obj = read(fidr, "obj_ref")
    @test obj.x.x === obj.x.y == obj.y.x === obj.y.y
    @test obj.x !== obj.y

    @check fidr padding_test
    @check fidr empty_arr_1
    @check fidr empty_arr_2
    @check fidr empty_arr_3
    @check fidr empty_arr_4
    @check fidr bigdata
    @check fidr bigfloats
    @check fidr bigints
    @check fidr none
    @check fidr nonearr
    @check fidr scalar_nothing
    @check fidr vector_nothing
    @check fidr Abig
    @check fidr Bbig
    @check fidr Sbig
    @check fidr bitsparamfloat
    @check fidr bitsparambool
    @check fidr bitsparamsymbol
    @check fidr bitsparamint
    @check fidr bitsparamuint
    @check fidr tuple_of_tuples

    x1 = read(fidr, "group1/x")
    @assert x1 == Any[1]
    x2 = read(fidr, "group2/x")
    @assert x2 == Any[2]

    close(fidr)
end
