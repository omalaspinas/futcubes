-- Inspired by the library of Paul Bourke

module mtetras = {

let tol = 1e-5f64

type point = { x:f64, y:f64, z:f64 }

let new_point (x': f64) (y': f64) (z': f64): point = {x = x', y = y', z = z'}

let zero: point = { x = 0.0f64, y = 0.0f64, z = 0.0f64 }

let sub (lhs: point) (rhs: point) : point = {
    x = lhs.x - rhs.x, 
    y = lhs.y - rhs.y, 
    z = lhs.z - rhs.z 
}

let cross (lhs: point) (rhs: point) : point = {
    x = lhs.y * rhs.z - lhs.z * rhs.y,
    y = - lhs.x * rhs.z + lhs.z * rhs.x,
    z = lhs.x * rhs.y - lhs.y * rhs.x
}

let norm_sqr (p: point) : f64 = 
    p.x * p.x + p.y * p.y + p.z * p.z

let norm (p: point) : f64 =
    f64.sqrt (norm_sqr p)

type triangle = {p0: point, p1: point, p2: point}

let empty_triangle : triangle = {
    p0 = zero,
    p1 = zero,
    p2 = zero
}

let area (t: triangle) : f64 = 
    let p01 = sub (t.p1) (t.p0)
    let p02 = sub (t.p2) (t.p0)

    in 0.5 * norm (cross p01 p02)

type grid_cell = { p: [8]point, value: [8]f64 }
 

let lookup (i: i64) : (i64, i64) = 
    -- Find the vertices where the surface intersects the cube 
    if i == 1 then (0, 1)
    else if i == 2 then (1, 2)
    else if i == 4 then (2, 3)
    else if i == 8 then (3, 0)
    else if i == 16 then (4, 5)
    else if i == 32 then (5, 6)
    else if i == 64 then (6, 7)
    else if i == 128 then (7, 4)
    else if i == 256 then (0, 4)
    else if i == 512 then (1, 5)
    else if i == 1024 then (2, 6)
    else if i == 2048 then (3, 7)
    else (-1, -1) -- use assert maybe

let vertex_interp (isolevel: f64) (p1: point) (p2: point) (valp1: f64) (valp2: f64): point =
    if f64.abs(isolevel - valp1) < tol then p1
    else if f64.abs(isolevel - valp2) < tol then p2
    else if f64.abs(valp1 - valp2) < tol then p1
    else
        let mu = (isolevel - valp1) / (valp2 - valp1)
        let px = p1.x + mu * (p2.x - p1.x)
        let py = p1.y + mu * (p2.y - p1.y)
        let pz = p1.z + mu * (p2.z - p1.z)
        in {x = px, y = py, z = pz}


-- Polygonise a tetrahedron given its vertices within a cube
-- This is an alternative algorithm to polygonisegrid.
-- It results in a smoother surface but more triangular facets.

--                  + 0
--                  /|\
--                 / | \
--                /  |  \
--               /   |   \
--              /    |    \
--             /     |     \
--            +-------------+ 1
--           3 \     |     /
--              \    |    /
--               \   |   /
--                \  |  /
--                 \ | /
--                  \|/
--                  + 2

-- It's main purpose is still to polygonise a gridded dataset and
-- would normally be called 6 times, one for each tetrahedron making
-- up the grid cell.
-- Given the grid labelling as in PolygniseGrid one would call
--     polygonise_tri (grid) (isolevel) (0) (2) (3) (7)
--     polygonise_tri (grid) (isolevel) (0) (2) (6) (7)
--     polygonise_tri (grid) (isolevel) (0) (4) (6) (7)
--     polygonise_tri (grid) (isolevel) (0) (6) (1) (2)
--     polygonise_tri (grid) (isolevel) (0) (6) (1) (4)
--     polygonise_tri (grid) (isolevel) (5) (6) (1) (4)
let polygonise_tri (grid: grid_cell) (isolevel: f64) (v0: i64) (v1: i64) (v2: i64) (v3: i64) : [2]triangle =
    -- Determine which of the 16 cases we have given which vertices
    -- are above or below the isosurface
    let tri_index = 0
    let tri_index = if grid.value[v0] < isolevel then tri_index | 1 else tri_index
    let tri_index = if grid.value[v1] < isolevel then tri_index | 2 else tri_index
    let tri_index = if grid.value[v2] < isolevel then tri_index | 4 else tri_index
    let tri_index = if grid.value[v3] < isolevel then tri_index | 8 else tri_index

    in 
        if tri_index == 0x0E || tri_index == 0x01 then 
            [
                { 
                    p0 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v1]) (grid.value[v0]) (grid.value[v1]),
                    p1 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v2]) (grid.value[v0]) (grid.value[v2]),
                    p2 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v3]) (grid.value[v0]) (grid.value[v3]) 
                }, 
                empty_triangle
            ]
            
        else if tri_index == 0x0D || tri_index == 0x02 then 
            [
                {
                    p0 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v0]) (grid.value[v1]) (grid.value[v0]),
                    p1 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v3]) (grid.value[v1]) (grid.value[v3]),
                    p2 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v2]) (grid.value[v1]) (grid.value[v2])
                },
                empty_triangle
            ]
        else if tri_index == 0x0C || tri_index == 0x03 then 
            let p1' = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v2]) (grid.value[v0]) (grid.value[v2])
            let p2' = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v3]) (grid.value[v1]) (grid.value[v3])
            in 
                [
                    {
                        p0 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v3]) (grid.value[v0]) (grid.value[v3]),
                        p1 = p1',
                        p2 = p2'
                    },
                    {
                        p0 = p2',
                        p1 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v2]) (grid.value[v1]) (grid.value[v2]),
                        p2 = p1'
                    }
                ]
        else if tri_index == 0x0B || tri_index == 0x04 then 
            [
                {
                    p0 = vertex_interp (isolevel) (grid.p[v2]) (grid.p[v0]) (grid.value[v2]) (grid.value[v0]),
                    p1 = vertex_interp (isolevel) (grid.p[v2]) (grid.p[v1]) (grid.value[v2]) (grid.value[v1]),
                    p2 = vertex_interp (isolevel) (grid.p[v2]) (grid.p[v3]) (grid.value[v2]) (grid.value[v3])
                },
                empty_triangle
            ]
        else if tri_index == 0x0A ||tri_index == 0x05 then
            let p0' = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v1]) (grid.value[v0]) (grid.value[v1])
            let p1' = vertex_interp (isolevel) (grid.p[v2]) (grid.p[v3]) (grid.value[v2]) (grid.value[v3])
            in 
                [
                    {
                        p0 = p0',
                        p1 = p1',
                        p2 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v3]) (grid.value[v0]) (grid.value[v3])
                    },
                    {
                        p0 = p0',
                        p1 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v2]) (grid.value[v1]) (grid.value[v2]),
                        p2 = p1'
                    }
                ]
        else if tri_index == 0x09 || tri_index == 0x06 then 
            let p0' = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v1]) (grid.value[v0]) (grid.value[v1])
            let p2' = vertex_interp (isolevel) (grid.p[v2]) (grid.p[v3]) (grid.value[v2]) (grid.value[v3])
            in 
                [
                    {
                        p0 = p0',
                        p1 = vertex_interp (isolevel) (grid.p[v1]) (grid.p[v3]) (grid.value[v1]) (grid.value[v3]),
                        p2 = p2'
                    },
                    {
                        p0 = p0',
                        p1 = vertex_interp (isolevel) (grid.p[v0]) (grid.p[v2]) (grid.value[v0]) (grid.value[v2]),
                        p2 = p2'
                    }
                ]
        else if tri_index == 0x07 || tri_index == 0x08 then 
            [
                {
                    p0 = vertex_interp (isolevel) (grid.p[v3]) (grid.p[v0]) (grid.value[v3]) (grid.value[v0]),
                    p1 = vertex_interp (isolevel) (grid.p[v3]) (grid.p[v2]) (grid.value[v3]) (grid.value[v2]),
                    p2 = vertex_interp (isolevel) (grid.p[v3]) (grid.p[v1]) (grid.value[v3]) (grid.value[v1])
                },
                empty_triangle
            ]
        else 
            [empty_triangle, empty_triangle]

let polygonise (grid: grid_cell) (isolevel: f64) : [12]triangle =
    flatten [
        (polygonise_tri (grid) (isolevel) (0) (2) (3) (7)) :> [2]triangle,
        (polygonise_tri (grid) (isolevel) (0) (2) (6) (7)) :> [2]triangle,
        (polygonise_tri (grid) (isolevel) (0) (4) (6) (7)) :> [2]triangle,
        (polygonise_tri (grid) (isolevel) (0) (6) (1) (2)) :> [2]triangle,
        (polygonise_tri (grid) (isolevel) (0) (6) (1) (4)) :> [2]triangle,
        (polygonise_tri (grid) (isolevel) (5) (6) (1) (4)) :> [2]triangle
    ] :> [12]triangle


let triangle_to_array (t: triangle) = 
    [
        t.p0.x, t.p0.y, t.p0.z,
        t.p1.x, t.p1.y, t.p1.z,
        t.p2.x, t.p2.y, t.p2.z
    ]


let polygonise_field [nx][ny][nz] (rho: [nx][ny][nz]f64) (isovalue: f64): []triangle = 
    let triangles = flatten_4d (
            tabulate_3d (nx-1) (ny-1) (nz-1) (\(x) (y) (z) ->
                let ooo: point = (new_point (f64.i64 (x))     (f64.i64 (y))     (f64.i64 (z)))
                let ioo: point = (new_point (f64.i64 (x + 1)) (f64.i64 (y))     (f64.i64 (z)))
                let iio: point = (new_point (f64.i64 (x + 1)) (f64.i64 (y + 1)) (f64.i64 (z)))
                let oio: point = (new_point (f64.i64 (x))     (f64.i64 (y + 1)) (f64.i64 (z)))

                let ooi: point = (new_point (f64.i64 (x))     (f64.i64 (y))     (f64.i64 (z + 1)))
                let ioi: point = (new_point (f64.i64 (x + 1)) (f64.i64 (y))     (f64.i64 (z + 1)))
                let iii: point = (new_point (f64.i64 (x + 1)) (f64.i64 (y + 1)) (f64.i64 (z + 1)))
                let oii: point = (new_point (f64.i64 (x))     (f64.i64 (y + 1)) (f64.i64 (z + 1)))

                let p': [8]point = [ooo, ioo, iio, oio, ooi, ioi, iii, oii]

                let value': [8]f64 = 
                    [
                        rho[i64.f64 (ooo.x), i64.f64 (ooo.y), i64.f64 (ooo.z)],
                        rho[i64.f64 (ioo.x), i64.f64 (ioo.y), i64.f64 (ioo.z)],
                        rho[i64.f64 (iio.x), i64.f64 (iio.y), i64.f64 (iio.z)],
                        rho[i64.f64 (oio.x), i64.f64 (oio.y), i64.f64 (oio.z)],
                        rho[i64.f64 (ooi.x), i64.f64 (ooi.y), i64.f64 (ooi.z)],
                        rho[i64.f64 (ioi.x), i64.f64 (ioi.y), i64.f64 (ioi.z)],
                        rho[i64.f64 (iii.x), i64.f64 (iii.y), i64.f64 (iii.z)],
                        rho[i64.f64 (oii.x), i64.f64 (oii.y), i64.f64 (oii.z)]
                    ]
                let grid: grid_cell = {p = p', value = value'}
                in (polygonise grid isovalue) :> [12]triangle
            )
        ) :> []triangle
    in filter (\t -> 
        area (t) > tol
    ) triangles -- only keep legal triangles
    
} -- module mcubes


