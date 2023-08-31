struct Element{T<:Real}
    a::T
    b::T
    c::T
    t::T
    x::T
    y::T
    izone::T
    id::T
    iphi::T
    sn::T
    co::T
    minr::T
    maxr::T
    minz::T
    maxz::T
end

Element(elements::AbstractMatrix{<:Real}, i::Int) = elm_data(@view(elements[:, i]))

function Element(i_data::AbstractVector{<:Real})
    a, b, c, t, x, y, _, izone, id, iphi = @views i_data[1:10]
    sn, co = sincos(t)
    ab = a + b
    x2 = x + ab * co
    x3 = x + b * co - c * sn
    y2 = y + ab * sn
    y3 = y + b * sn + c * co
    minr = min(x, x2, x3)
    maxr = max(x, x2, x3)
    minz = min(y, y2, y3)
    maxz = max(y, y2, y3)
    return Element(a, b, c, t, x, y, izone, id, iphi, sn, co, minr, maxr, minz, maxz)
end


function read_mesh(tid::HDF5.Group)
    mid = tid["mesh"]
    elements = read(mid, "elements")
    mesh = [Element(elm) for elm in eachcol(elements)]
    Nplanes = (read_attribute(mid, "3D") == 1) ? read_attribute(mid, "nplanes") : 0
    return mesh, Nplanes
end