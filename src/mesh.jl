struct Element{T<:Real}
    a::T
    b::T
    c::T
    t::T
    x::T
    y::T
    x2::T
    y2::T
    x3::T
    y3::T
    izone::T
    id::T
    iφ::T
    sn::T
    co::T
    minr::T
    maxr::T
    minz::T
    maxz::T
end

function localpos(r, z, localφ, elm)
    dr = r - elm.x
    dz = z - elm.y
    return @SVector[ dr * elm.co + dz * elm.sn - elm.b,
                    -dr * elm.sn + dz * elm.co,
                    localφ]
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



mutable struct BBnode{T<:Real}
    index::Int
    parent::Int
    left::Int
    right::Int
    minr::T
    maxr::T
    minz::T
    maxz::T
    area::T
    element::Union{Nothing, Element{T}}
end

function BBnode(tree::Vector{<:BBnode}, left, right; index=1, parent=0)
    # positive indices are
    minr = min(tree[left].minr, tree[right].minr)
    maxr = max(tree[left].maxr, tree[right].maxr)
    minz = min(tree[left].minz, tree[right].minz)
    maxz = max(tree[left].maxz, tree[right].maxz)
    area = (maxr - minr) * (maxz - minz)
    tree[left].parent  = index
    tree[right].parent = index
    return BBnode(index, parent, left, right, minr, maxr, minz, maxz, area, nothing)
end

function BBnode(elm::Element; index=1, parent=0)
    return BBnode(index, parent, 0, 0, elm.minr, elm.maxr, elm.minz, elm.maxz, area(elm), elm)
end

function update_extrema!(tree, node)
    BB = tree[node]
    left = tree[BB.left]
    right = tree[BB.right]
    BB.minr = min(left.minr, right.minr)
    BB.maxr = max(left.maxr, right.maxr)
    BB.minz = min(left.minz, right.minz)
    BB.maxz = max(left.maxz, right.maxz)
    BB.area = area(BB)
end

# function collision(base, test)
#     if ((base.minr <= test.minr <= base.maxr) && (base.minr <= test.maxr <= base.maxr) &&
#         (base.minz <= test.minz <= base.maxz) && (base.minz <= test.maxz <= base.maxz))
#         return :inside
#     elseif ((test.minr <= base.minr <= test.maxr) && (test.minr <= base.maxr <= test.maxr) &&
#             (test.minz <= base.minz <= test.maxz) && (test.minz <= base.maxz <= test.maxz))
#         return :outside
#     elseif (base.maxr >= test.minr) && (base.minr <= test.maxr) && (base.maxz >= test.minz) && (base.minz <= test.maxz)
#         return :overlap
#     else
#         return :none
#     end
# end
collision(base,test) = (base.maxr >= test.minr) && (base.minr <= test.maxr) && (base.maxz >= test.minz) && (base.minz <= test.maxz)

area(BB) = (BB.maxr - BB.minr) * (BB.maxz - BB.minz)
area(BB1::Union{BBnode, Element}, BB2::Union{BBnode, Element}) = (max(BB1.maxr, BB2.maxr) - min(BB1.minr, BB2.minr)) * (max(BB1.maxz, BB2.maxz) - min(BB1.minz, BB2.minz))

area(tree::Vector{<:BBnode}, node::Int) = tree[node].area
Δarea(tree::Vector{<:BBnode}, node::Int, leaf::Int) = area(tree[node], tree[leaf]) - area(tree, node)

function cost(tree, node, leaf)
    c = area(tree[node], tree[leaf])
    parent = tree[node].parent
    while parent != 0
        c += Δarea(tree, parent, leaf)
        parent = tree[parent].parent
    end
    return c
end

function bound_cost(tree, node, leaf)
    bc = area(tree, leaf)
    parent = node
    while parent != 0
        bc == Δarea(tree, parent, leaf)
        parent = tree[parent].parent
    end
    return bc
end

function find_sibling(tree, node, leaf)
    best = node
    mincost = cost(tree, node, leaf)
    if (tree[node].element === nothing) && (bound_cost(tree, node, leaf) < mincost)
        lbest, lcost = find_sibling(tree, tree[node].left,  leaf)
        rbest, rcost = find_sibling(tree, tree[node].right, leaf)
        if (lcost < mincost) && (lcost < rcost)
            best = lbest
            mincost = lcost
        elseif rcost < mincost
            best = rbest
            mincost = rcost
        end
    end
    return best, mincost
end

function get_sibling(tree, node)
    parent = tree[node].parent
    (parent == 0) && return nothing
    if tree[parent].left == node
        return tree[parent].right
    elseif tree[parent].right == node
        return tree[parent].left
    else
        throw(ErrorException("Node $(sibling) has parent $(parent), but node $(parent) does not have child $(sibling)"))
    end
end

function set_sibling!(tree, node, sibling)
    parent = tree[node].parent
    (parent == 0) && return nothing
    if tree[parent].left == node
        tree[parent].right = sibling
    elseif tree[parent].right == node
        tree[parent].left = sibling
    else
        println(tree)
        throw(ErrorException("Node $(node) has parent $(parent), but node $(parent) does not have child $(node)"))
    end
    tree[sibling].parent = parent
end

function rotate!(tree, node)
    parent = tree[node].parent
    parent == 0 && return
    grandpa = tree[parent].parent
    grandpa == 0 && return

    uncle   = get_sibling(tree, parent)
    sibling = get_sibling(tree, node)

    area1 = area(tree[node], tree[uncle])
    area2 = area(tree[sibling], tree[uncle])
    #println((area2, tree[parent].area))
    if (area2 < tree[parent].area) && (area2 <= area1)
        # Rotate
        #println("Rotating!")
        set_sibling!(tree, sibling, uncle)
        set_sibling!(tree, parent, node)
        update_extrema!(tree, parent)
    elseif area1 < tree[parent].area
        set_sibling!(tree, node, uncle)
        set_sibling!(tree, parent, sibling)
        update_extrema!(tree, parent)
    end
    rotate!(tree, parent)
end




function insert(tree, node, sibling, leaf)
    parent = tree[sibling].parent
    # create new node
    tree[node] = BBnode(tree, sibling, leaf; index=node, parent)

    # associate parent for the children
    tree[leaf].parent = node
    tree[sibling].parent = node

    # update child of the parent
    if parent != 0
        if tree[parent].left == sibling
            tree[parent].left = node
        elseif tree[parent].right == sibling
            tree[parent].right = node
        else
            throw(ErrorException("Node $(sibling) has parent $(parent), but node $(parent) does not have child $(sibling)"))
        end
    end

    while parent != 0
        update_extrema!(tree, parent)
        parent = tree[parent].parent
    end

    rotate!(tree, node)

end


struct BBtree{T}
    root::Int
    tree::Vector{BBnode{T}}
    mesh::Vector{Element{T}}
    Nplanes::Int
end

function BBtree(mesh::Vector{Element{T}}, Nplanes=0) where {T <: Real}
    N = length(mesh)
    (Nplanes > 0) && (N = N ÷ Nplanes)

    tree = Vector{BBnode{T}}(undef, 2N - 1)
    for k in 1:N
        tree[N+k-1] = BBnode(mesh[k]; index=k)
    end
    tree[1] = BBnode(tree, N, N+1; index=1, parent=0)
    root = 1
    for k in 2:(N-1)
        leaf = N + k
        best, _ = find_sibling(tree, root, leaf)
        insert(tree, k, best, leaf)
        (best == root) && (root = k)
    end
    return BBtree(root, tree, mesh, Int(Nplanes))
end

contains(node::BBnode, r, z) = (node.minr <= r <= node.maxr) && (node.minz <= z <= node.maxz)

function find_element(meshtree::BBtree, r::Real, z::Real, φ::Real=0.0, Nplanes::Integer=0)

    # get on proper plane
    offset = 0
    localφ = φ
    if Nplanes > 0
        N = length(meshtree.mesh)
        Np = N ÷ Nplanes
        while offset <= N
            elm = meshtree.mesh[offset+1]
            localφ = φ - elm.iφ
            if 0 <= localφ <= elm.id
                break
            else
                # this assumes elements are ordered by plane
                offset += Np
            end
        end
    end

    k = find_element(meshtree, meshtree.root, r, z, localφ)
    return (k == 0 ? 0 : k + offset), localφ
end

function find_element(meshtree::BBtree, k::Integer, r::Real, z::Real, localφ::Real=0.0)
    node = meshtree.tree[k]
    if contains(node, r, z)
        if node.element !== nothing
            lp = localpos(r, z, localφ, node.element)
            is_in_tri(lp, node.element) && return node.index
        else
            kl = find_element(meshtree, node.left, r, z, localφ)
            (kl > 0) && return kl
            kr = find_element(meshtree, node.right, r, z, localφ)
            (kr > 0) && return kr
        end
    end
    return 0
end

@recipe function plot_BBtree(tree::BBtree)
    tree, tree.root
end
@recipe function plot_BBtree(tree::BBtree, index)
    BB = tree.tree[index]
    @series begin
        linestyle --> :dash
        label := nothing
        color := :black
        @SVector[BB.minr, BB.maxr, BB.maxr, BB.minr, BB.minr], @SVector[BB.minz, BB.minz, BB.maxz, BB.maxz, BB.minz]
    end
    if BB.element === nothing
        @series begin
            label := nothing
            color := :red
            tree, BB.left
        end
        @series begin
            label := nothing
            color := :blue
            tree, BB.right
        end
    else
        BB.element
    end
end



Element(elements::AbstractMatrix{<:Real}, index::Int) = elm_data(@view(elements[:, i]); index)

function Element(i_data::AbstractVector{<:Real}; index = 0)
    a, b, c, t, x, y, _, izone, id, iφ = @views i_data[1:10]
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
    return Element(a, b, c, t, x, y, x2, y2, x3, y3, izone, id, iφ, sn, co, minr, maxr, minz, maxz)
end


function read_mesh(tid::HDF5.Group)
    mid = tid["mesh"]
    elements = read(mid, "elements")
    mesh = [Element(elm) for elm in eachcol(elements)]
    Nplanes = (read_attribute(mid, "3D") == 1) ? read_attribute(mid, "nplanes") : 0
    return mesh, Nplanes
end

@recipe function plot_element(elm::Element)
    @series begin
        seriestype --> :path
        linewidth --> 3
        @SVector[elm.x, elm.x2, elm.x3, elm.x], @SVector[elm.y,elm.y2, elm.y3, elm.y]
    end
end


@recipe function plot_mesh(mesh::Vector{<:Element}; height=600, pad=0.05)
    maxr = -Inf
    minr = Inf
    maxz = -Inf
    minz = Inf
    for elm in mesh
        minr = min(minr, elm.minr)
        maxr = max(maxr, elm.maxr)
        minz = min(minz, elm.minz)
        maxz = max(maxz, elm.maxz)
        (elm.iφ != 0.0) && break
        @series begin
            seriestype --> :path
            linewidth --> 1
            label := nothing
            @SVector[elm.x, elm.x2, elm.x3, elm.x], @SVector[elm.y,elm.y2, elm.y3, elm.y]
        end
    end
    rpad = pad * (maxr - minr)
    zpad = pad * (maxz - minz)
    legend --> false
    xlims --> (minr - rpad, maxr + rpad)
    ylims --> (minz - zpad, maxz + zpad)
    aspect_ratio --> :equal
    width = Int(ceil(height * sqrt((maxr - minr) / (maxz - minz))))
    size --> (width, height)
    xlabel --> "R"
    ylabel --> "Z"
    ()
end