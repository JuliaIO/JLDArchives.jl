using HDF5, JLD, Test

# names is deprecated to keys v0.12 JLD
if !isdefined(JLD, :keys)
   keys = names  
end

# Define variables of different types
x = 3.7
A = reshape(collect(1:15), 3, 5)
str = "Hello"
stringsA = String["It", "was", "a", "dark", "and", "stormy", "night"]
stringsU = String["It", "was", "a", "dark", "and", "stormy", "night"]
empty_string = ""
empty_string_array = String[]
empty_array_of_strings = String[""]
tf = true
TF = A .> 10
B = [-1.5 sqrt(2) NaN 6;
     0.0  Inf eps() -Inf]
AB = Any[A, B]
t = (3, "cat")
c = Float32(3)+Float32(7)im
cint = 1+im  # issue 108
C = reinterpret(ComplexF64, vec(B))
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
d = Dict([(syms[1], "aardvark"), (syms[2], "banana")])
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
vv = Vector{Int}[[1,2,3]] # issue #123
typevar = Array{Int}[[1]]
typevar_lb = Vector{<:Integer}[[1]]
typevar_ub = (Vector{U} where U>:Int)[[1]]
typevar_lb_ub = (Vector{U} where Int<:U<:Real)[[1]]
# Immutable type:
rng = 1:5
# Type with a pointer field (#84)
objwithpointer = r"julia"
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

iseq(x,y) = isequal(x,y)
iseq(x::MyStruct, y::MyStruct) = (x.len == y.len && x.data == y.data)
iseq(c1::Array{Base.Sys.CPUinfo}, c2::Array{Base.Sys.CPUinfo}) = length(c1) == length(c2) && all([iseq(c1[i], c2[i]) for i = 1:length(c1)])
function iseq(c1::Base.Sys.CPUinfo, c2::Base.Sys.CPUinfo)
    for n in Base.Sys.CPUinfo.names
        if getfield(c1, n) != getfield(c2, n)
            return false
        end
    end
    true
end
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
                error("For ", $(string(sym)), ", read value does not agree with written value")
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
checkexpr(a, b) = @test a == b
function checkexpr(a::Expr, b::Expr)
    @test a.head == b.head
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
    @test i >= length(a.args) && j >= length(b.args)
end

for fn in ("v0.2.26.jld", "v0.2.28.jld")
    for mmap = (true, false)
        undefv = Array{Any}(undef, 1)
        undefm = Array{Any}(undef, 2, 2)
        ms_undef = MyStruct(0)

        fidr = jldopen(joinpath(splitdir(@__FILE__)[1], fn), "r"; mmaparrays=mmap)
        @check fidr x
        !mmap && @check fidr A # FIXME: fails on 0.7 due to a data alignment issue
        @check fidr str
        @check fidr stringsA
        @check fidr stringsU
        @check fidr empty_string
        @check fidr empty_string_array
        @check fidr empty_array_of_strings
        @check fidr tf
        !mmap && @check fidr TF # FIXME: fails on 0.7 due to a data alignment issue
        !mmap && @check fidr AB # FIXME: fails on 0.7 due to a data alignment issue
        @check fidr t
        @check fidr c
        @check fidr cint
        !mmap && @check fidr C # FIXME: fails on 0.7 due to a data alignment issue
        @check fidr emptyA
        @check fidr emptyB
        !mmap && @check fidr ms # FIXME: fails on 0.7 due to a data alignment issue
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
        !mmap && @check fidr β # FIXME: fails on 0.7 due to a data alignment issue
        !mmap && @check fidr vv # FIXME: fails on 0.7 due to a data alignment issue
        @check fidr rng
        !mmap && @check fidr typevar # FIXME: fails on 0.7 due to a data alignment issue
        !mmap && @check fidr typevar_lb # FIXME: fails on 0.7 due to a data alignment issue
        !mmap && @check fidr typevar_ub # FIXME: fails on 0.7 due to a data alignment issue
        !mmap && @check fidr typevar_lb_ub # FIXME: fails on 0.7 due to a data alignment issue

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

        @test !in("objwithpointer", keys(fidr))
        @check fidr bt
        @check fidr sa_asc
        @check fidr sa_utf8
        @check fidr arr_empty_tuple

        x1 = read(fidr, "group1/x")
        @test x1 == Any[1]
        x2 = read(fidr, "group2/x")
        @test x2 == Any[2]

        # load but don't check
        read(fidr, "cpus")

        # subarray (changed representation in julia 0.4)
        if !mmap # FIXME: fails on 0.7 due to a data alignment issue
            subarray_safe = JLD.JLD00.readsafely(fidr, "subarray")
            @test subarray_safe["indexes"] == (1:5,)
            @test subarray_safe["parent"] == [1:5;]
        end

        close(fidr)
    end
end
