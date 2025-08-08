# BMP reference

## Data type and constructors
```@docs
BMP
BMP(::Integer, ::Integer)
BMP(::Integer, order)
projbmp
```

## Properties
```@docs
length(bmp::BMP)
bonddims
volume(bmp::BMP)
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
evalfunc(bmp::BMP, x::AbstractArray)
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
