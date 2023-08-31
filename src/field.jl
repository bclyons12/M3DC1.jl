const m  = @SVector[0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 3, 2, 1, 0]
const n  = @SVector[0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 2, 3, 4, 5]
const m1 = max.(m .- 1, 0)
const m2 = max.(m .- 2, 0)
const m3 = max.(m .- 3, 0)
const n1 = max.(n .- 1, 0)
const n2 = max.(n .- 2, 0)
const n3 = max.(n .- 3, 0)

const mi  = m  .+ 1
const mi1 = m1 .+ 1
const mi2 = m2 .+ 1
const mi3 = m3 .+ 1
const ni  = n  .+ 1
const ni1 = n1 .+ 1
const ni2 = n2 .+ 1
const ni3 = n3 .+ 1

const tmp = zeros(20)
const val = zeros(20)

function localpos(r, z, localphi, elm)
    dr = r - elm.x
    dz = z - elm.y
    return @SVector[ dr * elm.co + dz * elm.sn - elm.b,
                    -dr * elm.sn + dz * elm.co,
                    localphi]
end

is_in_tri(lp, elm::Element) = is_in_tri(lp, elm.a, elm.b, elm.c)
function is_in_tri(lp, a, b, c)
    small = (a + b + c) * 1e-4
    (lp[2] < 0 - small) && return false
    (lp[2] > c + small) && return false
    x = 1.0 - lp[2] / c
    (lp[1] < -b * x - small) && return false
    (lp[1] >  a * x + small) && return false
    return true
end

function eval_field(name::String, filename::String, slice::Int, r, z, phi=0.0; op::Int=1)
    fid = h5open(filename, "r")
    val = eval_field(name, fid, slice, r, z, phi; op)
    close(fid)
    return val
end

function eval_field(name::String, fid::HDF5.File, slice::Int, r, z, phi=0.0; op::Int=1)
    tid = fid[(slice == -1) ? "equilibrium" : "time_" * lpad(slice,3,"0")]
    return eval_field(name, tid, r, z, phi; op)
end

function eval_field(name::String, tid::HDF5.Group, r, z, phi=0.0; op::Int=1)
    field = read(tid["fields"], name)
    mesh, Nplanes = read_mesh(tid)
    return eval_field(field, mesh, r, z, phi; op, Nplanes)
end

function eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    Mr, Mz = LazyGrids.ndgrid(Rs, Zs)
   return _eval_field(field, mesh, Mr, Mz, phi; op, Nplanes)
end

function eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, Rs::AbstractVector{<:Real}, z::Real, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    Zs = [z for _ in eachindex(Rs)]
    return _eval_field(field, mesh, Rs, Zs, phi; op, Nplanes)
end

function eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, r::Real, Zs::AbstractVector{<:Real}, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    Rs = [r for _ in eachindex(Zs)]
    return _eval_field(field, mesh, Rs, Zs, phi; op, Nplanes)
end
function eval_field2(field::Matrix{<:Real}, mesh::Vector{<:Element}, r::Real, Zs::AbstractVector{<:Real}, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    #Rs = [r for _ in eachindex(Zs)]
    return _eval_field.(Ref(field), Ref(mesh), r, Zs, phi; op, Nplanes)
end

function eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, r::Real, z::Real, Phis::AbstractVector{<:Real}; op::Int=1, Nplanes::Integer=0)
    Rs = @SVector[r]
    Zs = @SVector[z]
    return _eval_field.(Ref(field), Ref(mesh), Ref(Rs), Ref(Zs), Phis; op, Nplanes)
end

function eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, r::Real, z::Real, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    Rs = @SVector[r]
    Zs = @SVector[z]
    return _eval_field(field, mesh, Rs, Zs, phi; op, Nplanes)[1]
end

function _eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, Rs::AbstractArray, Zs::AbstractArray, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    @assert size(Rs) === size(Zs)
    Vs = zeros(size(Rs))
    i = 0
    N = length(mesh)
    while i < N
        i += 1 # here so continue is handled properly
        elm = mesh[i]

        if Nplanes > 0
            localphi = phi - elm.iphi
            if (localphi < 0) || (localphi > elm.id)
                # this assumes elements are ordered by plane
                i = i + N ÷ Nplanes - 1
                continue
            end
        else
            localphi = zero(phi)
        end

        for k in eachindex(Vs)
            (Vs[k] != 0.0) && continue
            @inbounds z = Zs[k]
            ((z > elm.maxz) || (z < elm.minz)) && continue
            @inbounds r = Rs[k]
            ((r > elm.maxr) ||  (r < elm.minr)) && continue
            lp = localpos(r, z, localphi, elm)
            if is_in_tri(lp, elm)
                Vs[k] = eval(field, lp, elm.t, i; Nplanes, elm.izone, op, elm.sn, elm.co)
            end

        end
    end
    return Vs
end

function _eval_field(field::Matrix{<:Real}, mesh::Vector{<:Element}, r::Real, z::Real, phi::Real=0.0; op::Int=1, Nplanes::Integer=0)
    val = zero(promote_type(typeof(r), typeof(z), typeof(phi)))
    i = 0
    N = length(mesh)
    while i < N
        i += 1 # here so continue is handled properly
        elm = mesh[i]

        if Nplanes > 0
            localphi = phi - elm.iphi
            if (localphi < 0) || (localphi > elm.id)
                # this assumes elements are ordered by plane
                i = i + N ÷ Nplanes - 1
                continue
            end
        else
            localphi = zero(phi)
        end

        ((z > elm.maxz) || (z < elm.minz)) && continue
        ((r > elm.maxr) || (r < elm.minr)) && continue

        lp = localpos(r, z, localphi, elm)
        if is_in_tri(lp, elm)
            val = eval(field, lp, elm.t, i; Nplanes, elm.izone, op, elm.sn, elm.co)
            break
        end
    end
    return val
end


function eval(field, lp, t, ielm; Nplanes = 0, izone=1, op=1, sn = sin(t), co = cos(t))

    lp1 = @SVector[1., lp[1], lp[1]^2, lp[1]^3, lp[1]^4, lp[1]^5]
    lp2 = @SVector[1., lp[2], lp[2]^2, lp[2]^3, lp[2]^4, lp[2]^5]
    op2 = (op - 1) ÷ 10
    op1 = op - 10 * op2

    if op1 == 1
        tmp .= @. @views lp1[mi] * lp2[ni]
    elseif op1 == 2
        tmp .= @. @views co * m * lp1[mi1] * lp2[ni] - sn * n * lp1[mi] * lp2[ni1]
    elseif op1 == 3
        tmp .= @. @views sn * m * lp1[mi1] * lp2[ni] + co * n * lp1[mi] * lp2[ni1]
    elseif op1 == 4
        tmp .= @. @views (co * co * m * m1 * lp1[mi2] * lp2[ni] +
                         sn * sn * n * n1 * lp1[mi] * lp2[ni2] -
                         2 * co * sn * n * m * lp1[mi1] * lp2[ni1])
    elseif op1 == 5
        cs = co * sn
        tmp .= @. @views (cs * m * m1 * lp1[mi2] * lp2[ni] -
                         cs * n * n1 * lp1[mi] * lp2[ni2] +
                         (co * co - sn * sn) * n * m * lp1[mi1] * lp2[ni1])
    elseif op1 == 6
        tmp .= @. @views (sn * sn * m * m1 * lp1[mi2] * lp2[ni] +
                         co * co * n * n1 * lp1[mi] * lp2[ni2] +
                         2 * co * sn * n * m * lp1[mi1] * lp2[ni1])
    elseif op1 == 7
        tmp .= @. @views m * m1 * lp1[mi2] * lp2[ni] + n * n1 * lp2[ni2] * lp1[mi]
    elseif op1 == 8
        mm1 = m .* m1
        nn1 = n .* n1
        tmp .= @. @views (co * mm1 * m2 * lp1[mi3] * lp2[ni] -
                         sn * nn1 * n2 * lp2[ni3] * lp1[mi] -
                         sn * n * mm1 * lp1[mi2] * lp2[ni1] +
                         co * m * nn1 * lp1[mi1] * lp2[ni2])
    elseif op1 == 9
        mm1 = m .* m1
        nn1 = n .* n1
        tmp .= @. @views (sn * mm1 * m2 * lp1[mi3] * lp2[ni] +
                         co * nn1 * n2 * lp2[ni3] * lp1[mi] +
                         co * mm1 * n * lp1[mi2] * lp2[ni1] +
                         sn * m * nn1 * lp1[mi1] * lp2[ni2])
    end

    is3D = (Nplanes > 0)
    if op2 == 0
        @views @. val = field[1:20, ielm] + is3D * (field[21:40, ielm] * lp[3] +
                                                    field[41:60, ielm] * lp[3] ^ 2 +
                                                    field[61:80, ielm] * lp[3] ^ 3)
        tot = sum(val[k] * tmp[k] for k in 1:20)
    elseif op2 == 1 && is3D
        @views @. val = field[21:40, ielm] +
                        field[41:60, ielm] * lp[3] * 2.0 +
                        field[61:80, ielm] * (lp[3] ^ 2) * 3.0
        tot = sum(val[k] * tmp[k] for k in 1:20)
    elseif op2 == 2 && is3D
        @views @. val = field[41:60, ielm] * 2.0 +
                        field[61:80, ielm] * lp[3] * 6.0
        tot = sum(val[k] * tmp[k] for k in 1:20)
    else
        tot = zero(eltype(tmp))
    end
    return tot
end