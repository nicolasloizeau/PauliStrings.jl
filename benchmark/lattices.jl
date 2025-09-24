function xedges(Nx::Int, Ny::Int; bc="periodic")
    @assert bc in ["closed","periodic"]
    edges = []
    for y in 0:Ny-1
        for x in 0:Nx-2
            i = x%Nx + y *Nx + 1
            j = i+1
            push!(edges, (i,j))
        end
    end
    if bc == "periodic"
        for y in 0:Ny-1
            i = y *Nx + 1
            j = i+Nx-1
            push!(edges, (i,j))
        end
    end
    return edges
end

function yedges(Nx::Int, Ny::Int; bc="periodic")
    @assert bc in ["closed","periodic"]
    edges = []
    for x in 0:Nx-1
        for y in 0:Ny-2
            i = x + y *Nx + 1
            j = i+Nx
            push!(edges, (i,j))
        end
    end
    if bc == "periodic"
        for x in 0:Nx-1
            i = x + 1
            j = i+(Ny-1)*Nx
            push!(edges, (i,j))
        end
    end
    return edges
end

function square_lattice(Nx::Int, Ny::Int; bc="periodic")
    @assert bc in ["periodic", "closed", "yperiodic", "xperiodic"]
    edges = []
    xbc = "closed"
    ybc = "closed"
    if bc == "periodic" || bc == "xperiodic"
        xbc = "periodic"
    end
    if bc == "periodic" || bc == "yperiodic"
        ybc = "periodic"
    end
    return vcat(xedges(Nx, Ny; bc=xbc), yedges(Nx, Ny; bc=ybc))
end

"""
    lattice_to_site(a::Int, b::Int, Nx::Int, Ny::Int)

a: starts at 0
b: starts at 0
Nx, Ny: dimensions of the lattice
returns the site index starting at 1
"""
function lattice_to_site(a::Int, b::Int, Nx::Int, Ny::Int)
    @assert a<Nx && b<Ny
    i = b*Nx+a
    return i+1
end



function diagonal_lattice(Nx::Int, Ny::Int ; bc="closed")
    @assert bc in ["closed"]
    edges = []
    for a in 0:Nx-2
        for b in 0:Ny-2
            i = lattice_to_site(a, b, Nx, Ny)
            j = lattice_to_site(a+1, b+1, Nx, Ny)
            push!(edges, (i,j))
        end
    end
    for a in 1:Nx-1
        for b in 0:Ny-2
            i = lattice_to_site(a, b, Nx, Ny)
            j = lattice_to_site(a-1, b+1, Nx, Ny)
            push!(edges, (i,j))
        end
    end
    return edges
end


function plaquette_lattice(Nx::Int, Ny::Int; bc="closed")
    @assert bc in ["closed"]
    edges = []
    for a in 0:Nx-2
        for b in 0:Ny-2
            i = lattice_to_site(a, b, Nx, Ny)
            j = lattice_to_site(a+1, b, Nx, Ny)
            k = lattice_to_site(a, b+1, Nx, Ny)
            l = lattice_to_site(a+1, b+1, Nx, Ny)
            push!(edges, (i,j,k,l))
        end
    end
    return edges
end
