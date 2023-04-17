# RAD2DMod - Geometry Construction


## Initializing the model

Initializing an empty model starts with:
```@docs
create_empty_model()
```

## Basic Elements
There are several different functions for basic 2D geometry construction.

```@docs
edge(p1::Point2D, p2::Point2D; seed::T = 10, dir::Symbol=:pos) where T<:Integer
```