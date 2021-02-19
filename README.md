# Futcubes

A marching cube implementation in Futhark which is basically a port of
Paul Bourke's implementation of [Polygonising a scalar field](http://paulbourke.net/geometry/polygonise/).

## Installation

```bash
$ futhark pkg add github.com/omalaspinas/futcubes
$ futhark pkg sync
```

## Usage example

```futhark
import "lib/github.com/omalaspinas/futcubes/mcubes"

-- Returns an array of 9 f64 number which are the three
-- vertices of each triangle.

let main [nx][ny][nz] (rho: [nx][ny][nz]f64) (isovalue: f64): [][9]f64 = 
    map(\t -> mcubes.triangle_to_array t) (mcubes.polygonise_field rho isovalue)
```

