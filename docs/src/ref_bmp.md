# BMP reference

## Data type and constructors
```@docs
BMP
BMP(::Integer, ::Integer)
BMP(::Integer, ::Vector{<:Integer})
projbmp
```

## Properties
```@docs
length(bmp::BMP)
bonddims
volume
```

## Cleaning
```@docs
clean1_lr
clean1_rl
clean1
```

## Synthesis
```@docs
apply
minapply
multiapply
layerapply
```

## Other BMP operations
```@docs
evalfunc
insert_var
erase_var
restrict
joinfuncs
```

## Variable ordering
```@docs
swap!
reorder!
brute_force!
exact_minimize!
sift!
```
