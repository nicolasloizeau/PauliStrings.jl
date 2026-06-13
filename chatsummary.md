Here is a summary of the session transcript, detailing the problem investigated, the findings, the developed design, and where the work was interrupted.
The Problem

In the translation-symmetric (TS) case, the existing Trotter time evolution was slower than RK4. It expanded the TS Hamiltonian into all

        
O(N)
O(N)

      

translated ordinary Pauli terms, resulting in hundreds of sequential single-Pauli gate updates and frequent truncations per step. By contrast, TS RK4 operates directly on the compressed TS representation and performs only a few coarse-grained TS commutators.
Key Investigations & Lessons Learned

    The Allocation Hypothesis (Disproven): The initial assumption was that one-time resum(H) allocation and state expansion drove the slowdown. A prototype utilizing lazy gate generation (TrotterGatesTS) showed that avoiding this allocation did not resolve the bottleneck.

    The Granularity Bottleneck: Profiling showed the real culprit was the work granularity. Applying

            
    O(N)
    O(N)

          

    sequential gates—each requiring a dictionary update and truncation—created immense overhead compared to RK4's few global operations.

    Orbit-Level Splitting: To exploit TS, the Hamiltonian must be split into its

            
    m
    m

          

    translation orbits (

            
    H=∑Haorb
    H=∑Haorb​

          

    ) rather than individual translated terms, reducing the Trotter step to

            
    O(m)
    O(m)

          

    orbit-level flows instead of

            
    O(Nm)
    O(Nm)

          

    gates.

Designing the TS-Native Orbit Update

The challenge in orbit-level Trotter is that

        
eτLHa
eτLHa​

      

is not a simple single-Pauli rotation. The team developed and verified a TS Pauli-orbit rotation rule:

    Pauli-Orbit Graph: Applying the orbit Liouvillian to a TS Pauli string generates transitions to a small, closed set of other TS Pauli strings (average degree is small,

            
    ∼4
    ∼4

          

    nodes).

    Component Exponentials: This closed set forms a connected graph component. The exact orbit evolution of these strings is equivalent to a dense matrix exponential on this small component. In small-system tests, this matched exact dense evolution to machine precision.

Managing Component Count & Overhead

The exact component exponential was mathematically correct but slow because active operators decompose into hundreds of tiny components, making individual matrix exponentials and bookkeeping expensive. Two strategies were designed to address this:

    Component-Weight Freezing: Component weights are highly heavy-tailed; only

            
    ∼30
    ∼30

          

    out of

            
    ∼700
    ∼700

          

    components carry

            
    99.99%
    99.99%

          

    of the operator's weight. Freezing (carrying forward unchanged) low-weight components and only exponentiating high-weight ones yielded a

            
    ∼5–10×
    ∼5–10×

          

    speedup over RK4 on larger systems, though accuracy was limited by the second-order Strang splitting.

    Component Batching via Signatures: A diagnostic check revealed that many components are structurally identical (isomorphic). Grouping components by their "exact matrix signature" allows them to share a single computed matrix exponential (

            
    exp(τA)
    exp(τA)

          

    ) and enables batched matrix-matrix updates instead of individual matrix-vector multiplies.

Interruption Point

The user approved implementing the exact matrix signature batching method on the clean ts-active-orbit-trotter branch to dramatically reduce the overhead of the exact component calculations. However, the session terminated because the API budget limit for the gateway was reached.