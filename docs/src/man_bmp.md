# Basic BMP functionality
Here we introduce BMPs and various methods used to manipulate them. A `BMP` is a
collection of row-switching matrices (of type `RowSwitchMatrix`), a terminal
vector that contains the possible outputs of the product, and the variable
ordering information.

## Construction
There are a few ways of constructing BMPs. The most basic way is to use the
inner constructor of the `BMP` type, which allows the user to directly specify
its fields.
```julia
BMP(M::Matrix{RowSwitchMatrix}, R::Vector{<:Integer}, order::Vector{<:Integer})
```
Here, `M` should have size `(n, 2)` where `n` is the number of input bits. The
first dimension of this array corresponds to the position on the product, while
the second one corresponds to the value of the input bit at this position. The
`order` array must be of length `n`, consistent with `M`, and specifies the
label of the variable in each position. `R` is the terminal vector that
multiplies to the product.

In most cases, you should not use the inner constructor to create BMPs. The
starting points for most cases (e.g. when building the BMP for a circuit) are
constant-value BMPs or projection BMPs. The constant-value BMP is just the BMP
of a function that's identically zero or identically one. It can be obtained
using
```julia
BMP(val::Integer, n::Integer)
BMP(val::Integer, order::Vector{<:Integer})
```
where `val` is the value of the function. The second argument is either `n`, the
number of input variables, or `order`, which is variable ordering (from which
the number of input variables is inferred). Note that the variable ordering does
not change the resulting BMP in this case, but it will matter in the later
stages of BMP synthesis, so it should be set from the start.

A projection function is one that takes multiple inputs and returns one of them
as input. The interface provided in the library is as follows:
```julia
projbmp(xi::Integer, n::Integer)
projbmp(xi::Integer, order::Vector{<:Integer})
```
`xi` is the label (not the position) of the variable that's being projected
onto. The second argument follows the same conventions as in the constant-value
case.

### Input variables
Positions on the product are numbered `1` to `n`, going from left to right. The
variable labels also go from `1` to `n`, meaning that the `order` field must
contain integers in this range as well. You must ensure that this constraint is
satisfied when constructing BMPs, as otherwise you'll run into errors.

In manipulating BMPs you should always pay attention to which functions take
variable labels as inputs and which take variable positions. The default
behavior (when no ordering is specified in construction) is to order the
variables according to their labels, i.e. variable `1` is at site `1`, variable
`2` is at site `2`, etc. It is possible to modify the ordering after
construction as well, as explained later here.

### Utility functions
The library provides a few utility functions to access BMP data. The first of
these is `length`. This function from Julia base is overwritten for the `BMP`
type to return the number of input variables.
```julia
length(bmp::BMP)
```

The complexity of a BMP is characterized through its bond dimensions. These can
be obtained via
```julia
bonddims(bmp::BMP)
```
which returns an array containing the number of rows of each matrix along the
product. Alternatively, if you only need the size of a matrix at a particular
site you can use
```julia
bonddims(bmp::BMP, i::Integer)
```
where `i` is the position. The maximum bond dimension can be obtained using
```julia
max_dim(bmp::BMP)
```
Perhaps the most relevant quantity is what we refer to as the 'volume' of the
BMP, which is the sum of all bond dimensions (i.e. the row counts of the
matrices) and the terminal vector length. This number quantifies the total
amount of space the BMP uses.
```julia
volume(bmp::BMP)
```

### Evaluation
One of the most basic BMP operations is evaluation. This is simply the
computation of the output of the function represented by the BMP for given input
values.
```julia
evalfunc(bmp::BMP, x)
```
`x` should be an array of size `(n,...)` such that input sets are grouped along
the first dimension. The return value is an array of the same number of
dimensions. Its size along the first dimension is the number of output bits,
while the other dimensions match the input array's in size. (As an example, the
input array for a ``\{0,1\}^6 \to \{0,1\}^4`` function could be an `(6,3,2)`
array, in which case the size of the output would be `(4,3,2)`.)

Note that the input bits are ordered according to their labels and not
positions. This ensures that BMPs that are generated the same way but with
different variable orderings (either specified at the beginning or imposed
later) evaluate the same arrays to the same outputs. As a simple example,
consider the following:
```julia
all_inputs = [
    0 0 0 0 1 1 1 1;
    0 0 1 1 0 0 1 1;
    0 1 0 1 0 1 0 1;
] # (3,8) array

bmp1 = projbmp(3, [1,2,3])
out1 = evalfunc(bmp1, all_inputs)
@show out1 # (1,8) array

bmp2 = projbmp(3, [3,1,2])
out2 = evalfunc(bmp2, all_inputs)
@show out2 # (1,8) array
```
Both output arrays are the same.

### Bare BMPs
As metioned above, the `BMP` type stores some extra information about the
variable ordering and the terminal vector. In cases where these need to be
handled separately, or perhaps multiple disjoint BMPs with the same variable
ordering and terminal vectors are being processed, you may wish to use a
`BareBMP` which is just an alias for a `Matrix{RowSwitchMatrix}`. (This is
indeed how the `Chip` type is implemented.) Most of the BMP functionality is
defined on `BareBMP`, although you do need to use extra arguments in a few
cases.

## Compression or "cleaning"
An important BMP operation that is used internally by nearly all other
operations is called CLEAN. This operation compresses the BMP to its smallest
size possible with row-switching matrices. This operation can be performed
left-to-right (LTR) or right-to-left, to obtain different types of reductions.
In order to invoke these manually, you can use
```julia
clean1_lr(bmp::BMP) # LTR sweep
clean1_rl(bmp::BMP) # RTL sweep
```
There is also a convenience function, `clean1(bmp::BMP)`, which performs a LTR
sweep followed by a RTL sweep. Note that the default variants of the BMP
operations perform the necessary compression automatically, so you do not need
to call these functions yourself. However, you may want to delay this for
performance reasons, as CLEAN is a relatively costly operation. "No-cleaning"
variants of most functions are provided for such cases. It is explained below
what version of CLEAN should be preferred after each operation.

## Synthesis
The synthesis operation is called APPLY. This operation allows one to obtain the
BMP of a function
```math
h(f_1(\vec{x}), \dotsc, f_k(\vec{x}))
```
given the BMPs of ``f_1(\vec{x}), \dotsc, f_k(\vec{x})``. This is arguably the
most fundamental BMP operation and the chief method of obtaining BMPs of
complicated functions.

The most basic version of this operation can be called in one of two different
ways.
```julia
apply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
apply(bmps::Vector{BMP}, htab::Vector{<:Integer})
```
The first one is (essentially) a specialization of the second one for the
``k=2`` case. The ``k`` input BMPs are either specified in the array `bmps` or
in the arguments `bmp1` and `bmp2`. The final argument `htab` is a
one-dimensional array of length ``2^k``, containing the output value for each
input word. More explicitly, the value
```math
h(x_1, x_2, \dotsc, x_{k-1}, x_k)
```
is the ``(W+1)``-th element of the array `htab`, where
```math
W = 2^{k-1} x_1 + 2^{k-2} x_2 + \dotsc + 2 x_{k-1} + x_k.
```
Another version of APPLY called `minapply` exists with the same interface.
```julia
minapply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
minapply(bmps::Vector{BMP}, htab::Vector{<:Integer})
```
The BMPs generated by `apply` and `minapply` are always the same. The difference
is the algorithm used by these functions: `apply` is what is referred to as the
direct-product method, and `minapply` is the direct-sum method. If the bond
dimensions of the input BMPs are small, `apply` may be preferable, but for
larger BMPs, `minapply` will be faster. This is because `apply` generates larger
matrices before compression, while `minapply` avoids this with some extra
processing.

The difference between the two methods is clearer with the no-cleaning version
of each.
```julia
apply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
apply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
minapply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
minapply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
```
Note that these functions are not part of the public interface of the library,
so you need to import them manually as in
```julia
using BinaryMatrixProducts: apply_noclean
```
The function arguments work exactly the same way as before. If you need to
manually trigger compression after using these functions, you must use `clean1`
on any BMP obtained via `apply_noclean` and `clean1_rl` on any BMP obtained via
`minapply_noclean`. The algorithm used by `minapply` performs LTR cleaning "on
the fly", therefore there is no need for it later.

### An example
At this point we can provide complete examples of the library usage. In the
first one we generate the BMP of the function
```math
f(x_1, x_2, x_3) = \bar{x}_1 + x_2 x_3.
```
Here is one way of doing this.
```julia
using BinaryMatrixProducts

x1 = projbmp(1, 3)
x2 = projbmp(2, 3)
x3 = projbmp(3, 3)
const1 = BMP(1, 3)

t1 = apply(const1, x1, [0, 1, 1, 0])
t2 = apply(x2, x3, [0, 0, 0, 1])
bmp1 = apply(t1, t2, [0, 1, 1, 1])
```
We used the XOR operation to obtain the inversion of ``x_1``. In this example,
the final BMP can also be obtained instead in a single step as
```julia
bmp2 = apply([x1, x2, x3], [1, 1, 1, 1, 0, 0, 0, 1])
```
You can verify the equivalence of the two BMPs by evaluating them for various
input values using `evalfunc`.
```julia
all_inputs = [
    0 0 0 0 1 1 1 1;
    0 0 1 1 0 0 1 1;
    0 1 0 1 0 1 0 1;
]
@show all(evalfunc(bmp1, all_inputs) .== evalfunc(bmp2, all_inputs))
```
In addition, you can replace all occurrences of `apply` with `minapply` and
obtain the same final BMP.

A slightly more involved example is the BMP of an adder. This can be constructed
with the following function.
```julia
function build_adder_bmp(n::Integer)
    outputs = Vector{BMP}(undef, n+1)
    carry = BMP(0, 2*n)
    for i=1:n
        x = projbmp(i, 2*n)
        y = projbmp(i+n, 2*n)
        outputs[n+2-i] = apply([x, y, carry], [0, 1, 1, 0, 1, 0, 0, 1])
        carry = apply([x, y, carry], [0, 0, 0, 1, 0, 1, 1, 1])
    end
    outputs[1] = carry
    return outputs
end
```
This function returns an array of BMPs for the output bits of an `n` bit adder,
ordered from the most significant to the least significant.

### Working with joint BMPs
In the adder example, we created a BMP for each of the output bits. It is much
more convenient in general to have all the outputs in a single BMP. We can
obtain such a BMP from an array of disjoint BMPs using `joinfuncs`.
```julia
joinfuncs(bmps::Vector{BMP})
```
The resulting BMP has as outputs all the outputs of the individual BMPs stacked
on top of each other, in the order they're given.

Generating the output BMPs separately and combining them later may be preferable
in some cases. (If this is the case and you're dealing with reversible circuits,
you should consider using the `Chip` type.) In others, you may need to work with
a joint BMP from the start. When you do this, you need to use a different
synthesis function, called `multiapply`.
```julia
multiapply(
    bmp::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
```
`tabs` contains an array of truth tables following the same convention as
before, while the corresponding element of `bits` specifies which of the output
bits of `bmp` enter as inputs to a new output function. The resulting BMP
contains as outputs only the functions you specify. (In other words, if you want
to keep output `i` of `bmp` for further use, `bits` and `tabs` should contain
the elements `[i]` and `[0, 1]`, respectively. Note that the position of the
output will now be the position of these elements in the respective arrays, not
`i` necessarily.)

`multiapply` works with individual bitlines. Another function, `layerapply`, has
the same signature, but is designed as convenience for reversible circuits.
```julia
layerapply(
    bmp::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
```
A reversible circuit consists of reversible gates that have the same number of
outputs as inputs. Accordingly, `bits` is a list of gate inputs as before, but
`tabs` is now a permutation containing the output word for each input word. The
order of the bitlines is preserved: Regardless of where a gate is specified in
`bits` and `tabs`, all of its output bitlines will be in the same positions in
the new BMP as in `bmp`. It is assumed that all bitlines are part of exactly one
gate, and that the elements of `tabs` are permutations. You may get unexpected
results if these constraints are not satisfied.

## Ordering optimization
The last point that needs to be addressed is the dependence of BMP size on
variable ordering. This is a problem inherited from BDDs. The effect of variable
ordering can be drastic: There are examples for which one variable ordering
leads to BMPs of exponential size, while another ordering keeps them
linear-sized. The adder discussed above is such an example. In creating a BMP
for it, we kept the bits of the two operands packed together, which results in a
large BMP. A much better approach is to pair the corresponding bits of the
operands.

The optimization of ordering is an NP-complete problem. There are exact
solutions of exponential complexity, as well as fast heuristics that work well
in practice. Before we discuss these, let's just take a look at how we can
manually change the variable ordering. The basic method here is a local swap,
implemented by `swap!`.
```julia
swap!(bmp::BMP, i::Integer)
```
This function swaps the input variables on sites `i` and `i+1`. This change is
reflected in the `order` and `position` fields of the BMP as well. Note that `i`
specifies a position and not a variable label; the variable on site `i` moves to
site `i+1` and vice-versa. It is not difficult to see that all variable
orderings can be realized using local swaps. There is a convenience function for
this purpose.
```julia
reorder!(bmp, pord::Vector{<:Integer})
```
`pord` is the new variable ordering.

!!! note
    Unlike most of the functions explained in the previous sections, variable
    ordering functions modify the input BMP, as opposed to creating a new one.

!!! note
    Changing the variable ordering does not change the order in which variables
    should be specified in `evalfunc`. Evaluating a BMP with the same input bit
    array before and after a `swap!`, `reorder!`, or one of the other methods
    explained below must yield exactly the same result.

For the purpose of optimizing the ordering we provide two methods. One of these
is an exact minimization algorithm adapted from [the paper
here](https://ieeexplore.ieee.org/document/1512370). The other is heuristic
method widely used for BDDs, [developed originally by Ruddell in
1993](https://ieeexplore.ieee.org/document/580029), called sifting. The exact
mimization method uses a version of the A* algorithm enhanced with
branch-and-bound techniques, and is invoked by
```julia
exact_minimize!(bmp::BMP)
```
As mentioned earlier, this is an algorithm with exponential time complexity, and
it should not be used for large BMPs. The sifting algorithm should be the first
thing you try in most cases:
```julia
sift!(bmp::BMP)
```
This algorithm goes through the list of input variables, and optimizes the
position of each one while keeping others fixed. When it goes through all the
variables the algorithm is done. If you want this process to be repeated
multiple times, you may specify the number of such iterations.
```julia
sift!(bmp::BMP, n_iters::Integer)
```
You can try these methods on the adder.
```julia
n = 8
adder1 = build_adder_bmp(n)
@time exact_minimize!(adder1)
@show volume(adder1)
adder2 = build_adder_bmp(n)
@time sift!(adder2)
@show volume(adder2)
```
If you run this code, you will see that sifting is orders of magnitude faster.
It cannot find the optimal variable ordering, but the BMP it produces is still
within twice the size of the minimum.
