import "mcubes"

-- Toy example which should create an horizontal plane
-- at height 0.5
let dummy_example : []mcubes.triangle = 
    let ooo: mcubes.point = { x = 0f64, y = 0f64, z = 0f64 }
    let ioo: mcubes.point = { x = 1f64, y = 0f64, z = 0f64 }
    let iio: mcubes.point = { x = 1f64, y = 1f64, z = 0f64 }
    let oio: mcubes.point = { x = 0f64, y = 1f64, z = 0f64 }
    let ooi: mcubes.point = { x = 0f64, y = 0f64, z = 1f64 }
    let ioi: mcubes.point = { x = 1f64, y = 0f64, z = 1f64 }
    let iii: mcubes.point = { x = 1f64, y = 1f64, z = 1f64 }
    let oii: mcubes.point = { x = 0f64, y = 1f64, z = 1f64 }

    let value': [8]f64 = [0f64, 0, 0, 0, 1, 1, 1, 1]
    let p': [8]mcubes.point = [ooo, ioo, iio, oio, ooi, ioi, iii, oii]

    let grid: mcubes.grid_cell = {p = p', value = value'}
    in filter (\t -> mcubes.area t > mcubes.tol) (mcubes.polygonise grid 0.5)

-- ==
-- entry: main
-- input { }
-- output { [[1.000000f64, 0.000000f64, 0.500000f64, 0.000000f64, 0.000000f64, 0.500000f64, 1.000000f64, 1.000000f64, 0.500000f64], 
--          [1.000000f64, 1.000000f64, 0.500000f64, 0.000000f64, 0.000000f64, 0.500000f64, 0.000000f64, 1.000000f64, 0.500000f64]] }
entry main : [][9]f64 =
    map (mcubes.triangle_to_array) (dummy_example)

-- ==
-- entry: main_field
-- input { [0f64, 1f64, 0f64, 1f64, 0f64, 1f64, 0f64, 1f64] }
-- output { [[1.000000f64, 0.000000f64, 0.500000f64, 0.000000f64, 0.000000f64, 0.500000f64, 1.000000f64, 1.000000f64, 0.500000f64], 
--          [1.000000f64, 1.000000f64, 0.500000f64, 0.000000f64, 0.000000f64, 0.500000f64, 0.000000f64, 1.000000f64, 0.500000f64]] }
entry main_field [n] (rho: [n]f64 ) : [][9]f64 =
    let rho' = unflatten_3d 2 2 2 rho
    in map (mcubes.triangle_to_array) (mcubes.polygonise_field rho' 0.5)