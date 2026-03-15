# Example demonstrating finding LIOMs in the XXZ model. Follows section III.A of https://arxiv.org/abs/2505.05882

using PauliStrings
import PauliStrings as ps
using Printf
using PyPlot

function XXZ(N, J, Δ)
    H = ps.Operator(N)
    i = 1
    H += J / 2, "S+", i, "S-", mod1(i + 1, N)
    H += J / 2, "S-", i, "S+", mod1(i + 1, N)
    H += J * Δ, "Sz", i, "Sz", mod1(i + 1, N)
    return ps.OperatorTS1D(H, full=false)
end

Δ = 0.75
J = 1.0

# ======================== Symmetry-adapted basis LIOMs ==========================

symmetry_Ms = [2, 4, 6]
symmetry_mazurs = []
symmetry_evals = []

for M in symmetry_Ms
    L = 2 * M + 1
    H = XXZ(L, J, Δ)

    support = ps.symmetry_adapted_k_local_basis(
        L,
        M;
        time_reversal=:real,
        spin_flip=:even,
        conserve_magnetization=:yes,
        translational_symmetry=true,
    )

    evals, evecs, ops = ps.lioms(
        H,
        support;
        threshold=Inf # return all eigenmodes
    )
    push!(symmetry_evals, evals)

    if M <= 4
        lioms = ops[evals.<1e-10]
        mazurs = zeros(Float64, length(support))
        for (i, liom) in enumerate(lioms)
            for (j, op) in enumerate(support)
                # mazurs[j] += (ps.trace_product(dagger(liom), op) / ps.trace_product(dagger(op), op))^2
                mazurs[j] += evecs[j, i]^2
            end
        end
        push!(symmetry_mazurs, mazurs)
    end
end


# ================================ Plotting ================================

fig, (ax1, ax2) = subplots(1, 2, figsize=(10, 4))

markers = ["o", "s", "^", "v", "D", "<", ">"]
for (i, M) in enumerate(symmetry_Ms)
    evals = sort(symmetry_evals[i])
    α = collect(1:length(evals))
    ax1.plot(α, evals, marker=markers[mod1(i, length(markers))],
        label="M = $M", markersize=5, linewidth=1)
end

ax1.set_xlabel("α", fontsize=12)
ax1.set_ylabel(L"Eigenvalue $\lambda_{\alpha}$", fontsize=12)
ax1.set_xlim(0.5, 9.5)
ax1.set_ylim(-0.02, 0.52)
ax1.set_xticks(1:1:9)
ax1.legend(loc="upper right", fontsize=11)
ax1.grid(true, alpha=0.3)

# Mazur
for (i, M) in enumerate(symmetry_Ms)
    if M > 4
        continue
    end
    mazur = symmetry_mazurs[i]
    s = collect(1:length(mazur))
    ax2.bar(s .+ (i - 2) * 0.5, mazur, width=0.4,
        label="M = $M", alpha=0.7, edgecolor="black", linewidth=0.8)
end
ax2.legend(loc="upper right", fontsize=12)
ax2.set_xlabel(L"$L$", fontsize=11)
ax2.set_ylabel(L"Mazur Bound $B_{s}$", fontsize=12)
ax2.set_xlim(0.5, 17)
ax2.set_ylim(1e-6, 1.2)
ax2.set_xticks(1:1:16)
ax2.set_yscale("log")

tight_layout()
savefig("symmetry_resolved_lioms_xxz.png", dpi=300)
