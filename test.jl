using Meshes, MeshViz, Rotations
import GLMakie
using LinearAlgebra

# verticies of a cube
v = Point3f[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

# faces - open on top and one side
f = [(1, 2, 3), (3, 4, 1), (1, 2, 5), (5, 6, 2), (2, 3, 7), (7, 6, 2), (3, 4, 8), (8, 7, 3)]


# mesh
c = connect.(f, Triangle)
m = SimpleMesh(v, c)


# save to stl
function writeSTL_ascii(mesh::Mesh, filename::String)
    file = filename * ".stl"
    N = length(mesh)
    facets = join(map(i -> write_facet_ascii(get_facet_data_from_mesh(mesh, i)...), 1:N), "\n\t")
    str = """
          solid name
              $facets
          endsolid name        
          """
    open(file, "w") do io
        write(io, str)
    end

end

# write an individual facet
function write_facet_ascii(n, v)
    str = """
          facet normal $(n[1]) $(n[2]) $(n[3])
              outer loop
                  vertex $(v[1].coords[1]) $(v[1].coords[2]) $(v[1].coords[3])
                  vertex $(v[2].coords[1]) $(v[2].coords[2]) $(v[2].coords[3])
                  vertex $(v[3].coords[1]) $(v[3].coords[2]) $(v[3].coords[3])
              endloop
          endfacet
          """

    return str
end

# get the indicies of the vertices of a facet
get_idxs_from_mesh(m::Mesh, idx) = m.topology.connec[idx].indices

# get the normal of a facet
get_normal_from_vs(v1, v2, v3) = normalize((v2 - v1) × (v3 - v1))

# get the normal and vertices of a facet
function get_facet_data_from_mesh(mesh::Mesh, idx)
    idx1, idx2, idx3 = get_idxs_from_mesh(mesh, idx)
    v = [mesh.vertices[idx1], mesh.vertices[idx2], mesh.vertices[idx3]]
    n = get_normal_from_vs(v...)

    return n, v
end


# write binary stl file
#=
specification:
    - 80 byte header, generally ignored but should not contain `solid` to be confused with ascii stl [80 UINT8]
    - 4 byte unsigned integer (little endian), number of facets [UINT32]
    - triangle data, each one described by 12 32-bit floating point numbers (little endian)
        - 3 32-bit floating point numbers, normal vector [3 REAL32]
        - 3 32-bit floating point numbers, vertex 1 [3 REAL32]
        - 3 32-bit floating point numbers, vertex 2 [3 REAL32]
        - 3 32-bit floating point numbers, vertex 3 [3 REAL32]
        - 2 byte unsigned integer, attribute byte count (always 0) [UINT16]
=#

function writeSTL_binary(mesh::Mesh, filename)
    file = filename * ".stl"

    open(file, "w") do io
        # header
        write(io, Vector{UInt8}("UNITS=mm"))
        write(io, zeros(UInt8, 72))

        # number of facets
        N = length(mesh) |> UInt32
        write(io, N)

        # triangle data
        for i in 1:N
            # normal
            n, v = get_facet_data_from_mesh(mesh, i)
            write(io, Float32.(n)...)

            # vertices
            foreach(x -> write(io, Float32.(x.coords)...), v)

            # attribute byte count
            write(io, zeros(UInt16, 1))
        end
    end

end



##### can save with MeshIO as STL but need to provide necessary conversions #####


### gull wing ###
# functions from above
camber(η, S, zcmax) = zcmax * η * (η - 1) * sum(S .* (2η - 1) .^ (0:2))
thickness(η, A, ztmax) = ztmax * sum(A .* (η .^ (2:5) .- √η))
chord(ζ, E, c₀,) = c₀ * (Fok(ζ) .+ Fcorr(ζ, E))
Fok(ζ) = 0.0 ≤ ζ ≤ 0.5 ? 1.0 : 4ζ * (1 - ζ)
Fcorr(ζ, E) = sum(@.E * (ζ^(3:7) - ζ^8))

S = [3.8735, -0.807, 0.771]
A = [-15.246, 26.482, -18.975, 4.6232]
E = [26.08, -209.92, 637.21, -945.21, 695.03]
c₀ = 0.388
zcmax(ζ) = 0.14 / (1 + 1.333ζ^1.4)
ztmax(ζ) = 0.1 / (1 + 3.546ζ^1.4)


function wing_mesh(S, A, zcmax, ztmax, E, c₀)
    # mesh parameters
    Nchord = 10 # number of points per aerofoil
    Nspan = 10 # number of aerofoils

    # spanwise and chordwise locations
    ζs = range(0, 1, length=Nspan)
    ηs = range(0, 1, length=Nchord ÷ 2)

    # xyz points 
    points = Meshes.Point3[]

    for ζ in ζs
        # camber line
        camber_line = [camber(η, S, -zcmax(ζ)) for η in ηs]

        # thickness distribution
        thickness_dist = [thickness(η, A, ztmax(ζ)) for η in ηs]

        # upper and lower surfaces of aerofoil
        upper_surface = camber_line + thickness_dist
        lower_surface = camber_line - thickness_dist
        # surface = [upper_surface; lower_surface]

        # aerofoils points at current spanwise coordinate
        # using aircraft coordinate system: x forward, y spanwise, z vertical
        aerofoil1 = map(zip(ηs, upper_surface)) do (η, p)
            Meshes.Point3(η * chord(ζ, E, c₀), ζ, p * chord(ζ, E, c₀))
            # __v = [η * chord(ζ, E, c₀), 0.0, p * chord(ζ, E, c₀)]
            # __v = ζ ≥ 0.0 ? Vector(RotX(deg2rad(-42)) * __v) : __v
            # __v[2] += ζ
            # Meshes.Point3(__v)

        end
        aerofoil2 = map(zip(ηs, lower_surface)) do (η, p)
            Meshes.Point3(η * chord(ζ, E, c₀), ζ, p * chord(ζ, E, c₀))
            # __v = [η * chord(ζ, E, c₀), 0.0, p * chord(ζ, E, c₀)]
            # __v = ζ ≥ 0.0 ? Vector(RotX(deg2rad(-42)) * __v) : __v
            # __v[2] += ζ
            # Meshes.Point3(__v)

        end

        aerofoil = [aerofoil1; reverse(aerofoil2)]


        append!(points, aerofoil)
    end

    # bend lower arm
    Npoints = length(points)
    points[(Npoints÷2):end] = map(points[(Npoints÷2):end]) do point
        rot = recenter(LinearMap(RotZ(-π / 4)), [0.0, 0.5, 0.0])
        temp = rot(point.coords) |> Vector

        return Point3f(temp...)
    end




    # connections
    faces = NTuple{3,Int64}[]
    for i in 1:Nspan-1
        for j in 1:Nchord-1
            a = (i - 1) * Nchord + j # cross-section 1, point 1
            b = a + 1 # cross-section 1, point 2
            c = a + Nchord # cross-section 2, point 1
            d = c + 1

            # push!(faces, (a, b, d, c))
            push!(faces, (a, b, c))
            push!(faces, (c, d, b))
        end
    end

    # triangulate the end face



    # end_faces = [tuple(1:Nchord...), tuple(1+Nspan:Nchord+Nspan...)]

    # connections = connect.([faces; end_faces], Ngon)
    connections = connect.(faces, Triangle)


    # create mesh
    return SimpleMesh(points, connections)
end
