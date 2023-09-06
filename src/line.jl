struct Line{T}
    r1::T
    φ1::T
    z1::T
    r2::T
    φ2::T
    z2::T
end

Line(r1, φ1, z1, r2, φ2, z2) = Line(promote(r1, φ1, z1, r2, φ2, z2)...)

function parametric_point(L::Line, t::Real)
    r = L.r1 + (L.r2 - L.r1) * t
    z = L.z1 + (L.z2 - L.z1) * t
    φ = L.φ1 + (L.φ2 - L.φ1) * t
    return r, z, φ
end

function Base.length(L::Line)
    x1 = L.r1 * cos(L.φ1)
    y1 = L.r1 * sin(L.φ1)
    x2 = L.r2 * cos(L.φ2)
    y2 = L.r2 * sin(L.φ2)
    return sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (L.z2 - L.z1) ^ 2)
end

function eval_field(name::String, tid::HDF5.Group, meshtree::BBtree, L::Line, Npts::Int; op::Int=1)
    field = read(tid["fields"], name)
    return eval_field(field, meshtree, L, Npts; op)
end

function eval_field(field::Matrix{<:Real}, meshtree::BBtree, L::Line, Npts::Int; op::Int=1)
    return [_eval_field(field, meshtree, parametric_point(L, t)...; op) for t in range(0, 1, Npts)]
end

function line_integral(name::String, tid::HDF5.Group, meshtree::BBtree, L::Line, Npts::Int; op::Int=1)
    field = read(tid["fields"], name)
    return line_integral(field, meshtree, L, Npts; op)
end

function line_integral(field::Matrix{<:Real}, meshtree::BBtree, L::Line, Npts::Int; op::Int=1)
    V = eval_field(field, meshtree, L, Npts; op)
    dl = length(L) / (Npts - 1)
    return @views dl * (0.5 *(V[1] + V[end]) + sum(V[2:(end-1)]))
end