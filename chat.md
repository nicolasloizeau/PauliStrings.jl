 ▐▛███▜▌   Claude Code v2.1.167
▝▜█████▛▘  gpt-5.5 · API Usage Billing
  ▘▘ ▝▝    ~/Programs/unitaryHack/PauliStrings.jl


❯ hi                                                                            

● Hi! What would you like to work on in PauliStrings.jl?

✻ Churned for 2s

❯ Context                                                                       
  The current Trotter time evolution function for translation symmetric         
  operators expands the Hamiltonian via resum(H), Trotterizes the full          
  operator, and evolves the representative of O through it. This is correct but 
  does not outperform RK4, defeating the point of Trotter, which beats RK4 in   
  the non-TS case.                                                              
  This is exemplified in https://github.com/nicolasloizeau/PauliStrings.jl/blob 
  /main/examples/evolve_trotter_ts.jl                                           
                                                                                
  Why does the current resum-then-Trotter approach lose advantage in the TS     
  case ? Design a second order Trotter method that exploits translation         
  invariance and outperform RK4 in the TS case.                                 
  The above example is a minimal fast running code.  For performance            
  benchmarks, you can use the same mixed field ising model at longer time and   
  larger system size.                                                           
                                                                                
  This is your problem. First find what is causing the slowdown, and why the    
  non-TS case uses the optimized codepath. Then, design. Be sure to cite all    
  external works. If necessary profile things, but keep in mind this can be     
  slow. Do not implement these yet. Let us first discuss the problem            

● Explore(Explore Trotter code) Sonnet 4.6
  ⎿  Done (0 tool uses · 0 tokens · 0s)
  (ctrl+o to expand)

  Searched for 3 patterns, read 12 files (ctrl+o to expand)

  Read 1 file, ran 5 shell commands (ctrl+o to expand)

● Web Search("Trotter 1959 product formula citation Proc Amer Math Soc On the
            product of semi-groups of operators")
  ⎿  Did 0 searches in 1s

❯ what are you up to                                                            

● I’m investigating the slowdown and sketching a TS-aware Trotter design,
  without implementing anything yet.

  Current findings:

  - The TS evolve(...; method=Trotter()) path is in src/evolution.jl:247-279.
  - It does:
    a. Or = representative(O)
    b. Hr = resum(H)
    c. gates = trotterize(Hr, dt)
    d. trotter_step!(Or, gates; truncation=...)
    e. wraps output back as OperatorTS{Ls,Ps}(Or) only for saving/final return.

  So the algorithm throws away most of the TS advantage during the inner loop.

  For the mixed-field Ising example:

  H = OperatorTS{(N,)}(3 local terms)
  Hr = resum(H)  # 3N full Pauli terms

  Then second-order Strang Trotter creates roughly:

  2 * length(Hr) - 1 = 6N - 1

  single-Pauli rotation gates per time step.

  So at N=64, the TS Hamiltonian has only 3 representative terms, but current
  Trotter applies 383 gates per step.

  By contrast, TS RK4 keeps using Operator{<:PauliStringTS} arithmetic. Its
  commutator kernel in src/translation_symmetry.jl:260-289 loops over
  representative terms and shifts, then recompresses into TS representatives.
  That means RK4 pays for translation sums inside the TS algebra but keeps the
  evolving operator compressed.

  Quick timing probes, not a full benchmark:

  - N=32, M=2^12, 50 steps: RK4 ≈ 1.37s, current TS Trotter ≈ 1.33s
  - N=64, M=2^14, 20 steps: RK4 ≈ 2.49s, current TS Trotter ≈ 4.02s
  - H length: 3 TS terms → resum(H) length 3N
  - Strang gates: 6N - 1

  The non-TS Trotter advantage comes from src/trotter.jl:75-109: each gate is a
  single Pauli rotation applied directly to each current Pauli string, using a
  cheap commutator and 2-branch update. That avoids the 4 full commutator
  evaluations and intermediate additions of RK4. But in the TS case, resum(H)
  multiplies the gate count by N, and evolving representative(O) as a plain
  Operator loses the compression that makes TS RK4 competitive.

  The design direction I’m considering:

  - Do not build trotterize(resum(H)).
  - Build Trotter layers from the representative local Hamiltonian terms.
  - Apply each translated local term to the TS operator while staying in
  OperatorTS.
  - Use a second-order Strang structure over translation-invariant Hamiltonian 
  components, e.g. for MFIM:

  exp(dt/2 L_Hx) exp(dt/2 L_Hz) exp(dt L_Hzz)
  exp(dt/2 L_Hz) exp(dt/2 L_Hx)

  - where each L_Ha represents the sum over all translations of one local term.

  The key unresolved design issue is how to apply one TS Trotter factor
  efficiently:

  exp(i dt ∑_r h_r) O exp(-i dt ∑_r h_r)

  without expanding into N independent gates in the hot loop. There are two
  plausible designs:

  1. Layered commuting-orbit rotations, when all translations in a Hamiltonian
  component commute.
  This is ideal for onsite X, onsite Z, and Ising ZZ terms: all translated terms
  within each component commute, so their product can be applied as a
  translation orbit. We can apply an orbit gate to OperatorTS and recompress
  after each orbit/layer.
  2. Even/odd or graph-colored local layers, for non-commuting translated copies
  of the same local term.
  For nearest-neighbor terms, split bonds into commuting color classes. In 1D
  nearest-neighbor periodic chains, usually two colors suffice for even/odd
  bonds, with care for odd N. Then Strang-symmetrize the color sequence.

  My tentative recommendation: start with a TS-specific second-order Trotter for
  Hamiltonians represented as sums of local TS terms, using orbit/layer gates.
  For the mixed-field Ising benchmark this should reduce per-step gates from O(N
  * n_terms) to roughly O(n_terms) TS orbit applications, while preserving TS
  compression of O. That is the path likely to beat RK4.

✻ Baked for 38s

❯ /resume 
  ⎿  Resume cancelled
※ recap: We’re analyzing why translation-symmetric Trotter is slow and designing a TS-aware second-order method before implementation. Next I’ll finish the design comparison and identify the most promising algorithm to discuss. (disable recaps in /config)

❯ Always try to take paths which are generic and do not rely on additional      
  structure of the problem                                                      

  Read 1 file (ctrl+o to expand)

● Fetch(https://arxiv.org/abs/2505.21606)
  ⎿  Received 44.6KB (200 OK)

● Fetch(https://epubs.siam.org/doi/10.1137/0705041)
  ⎿  Received 1.7MB (200 OK)

  Recalled 1 memory, wrote 2 memories (ctrl+o to expand)

● Noted — that’s the right constraint. I’ll make the main design generic, and
  treat things like “MFIM has commuting onsite/ZZ layers” or “1D
  nearest-neighbor even/odd coloring” only as optional accelerations, not as the
  core method. I also saved this as a working preference for future design
  discussions.

● Fetch(https://arxiv.org/abs/math-ph/0506007)
  ⎿  Received 45.3KB (200 OK)

● Fetch(https://www.ams.org/journals/proc/1959-010-04/S0002-9939-1959-0108732-6/
       S0002-9939-1959-0108732-6.pdf)

● Fetch(https://arxiv.org/pdf/2505.21606)

  Ran 1 shell command (ctrl+o to expand)

❯ /resume 
  ⎿  Resume cancelled

❯ continue                                                                      

● Agreed. With that constraint, I would not make the main solution depend on
  MFIM-specific commuting layers, even/odd bond coloring, nearest-neighbor
  geometry, etc. Those can be optional later optimizations. The generic path is:

  ▎ keep the mathematically standard full Strang product over translated Pauli 
  ▎ terms, but make the state representation and truncation 
  ▎ translation-symmetric after every complete Trotter step.

  That is generic for any OperatorTS Hamiltonian.

  What is causing the slowdown?

  Current TS Trotter path

  The translation-symmetric Trotter implementation is in
  src/evolution.jl:247-279.

  It currently does:

  Or = representative(O)
  Hr = resum(H)
  gates = trotterize(Hr, dt)
  trotter_step!(Or, gates; truncation=truncation)
  _save!(..., OperatorTS{Ls,Ps}(Or), ...)

  Relevant lines:

  - resum(H) expands the TS Hamiltonian into all translated ordinary strings:
  src/evolution.jl:260
  - gates are built from that expanded ordinary Hamiltonian:
  src/evolution.jl:265-273
  - trotter_step! evolves Or as an ordinary Operator: src/evolution.jl:274
  - the result is wrapped as OperatorTS only for saving/final return:
  src/evolution.jl:277,279

  resum itself is simple full unfolding:

  for s in all_shifts(paulistringtype(o))
      op += shift(rep_op, Ls, Ps, s)
  end

  at src/translation_symmetry.jl:248-258.

  For the example Hamiltonian:

  H = OperatorTS{(N,)}(3 representative local terms)

  resum(H) has 3N ordinary Pauli terms. Then second-order Strang gives about:

  2 * 3N - 1 = 6N - 1

  ordinary single-Pauli rotation gates per step.

  For N = 64, that is 383 gates per step from only 3 TS representative
  Hamiltonian terms.

  Why RK4 does better in the TS case

  RK4 is generic over AbstractOperator, so with H::OperatorTS and O::OperatorTS,
  the commutator uses the TS binary kernel in
  src/translation_symmetry.jl:260-289.

  That kernel does not expand both operands into full ordinary operators. It
  loops over:

  representative terms of A
  representative terms of B
  relative translations

  and immediately canonicalizes results back into PauliStringTS.

  So RK4 pays for translations inside the kernel, but the state stays compressed
  in the TS basis.

  The current Trotter path, by contrast, starts from a representative ordinary
  operator, applies all full translated gates, and then keeps evolving that
  ordinary operator. It does not re-compress to TS before the next step.

  This causes two related problems:

  1. Gate count inflation: H goes from m TS terms to mN ordinary terms.
  2. State inflation: O is evolved as an ordinary Operator, so
  translated-equivalent strings are not folded together between steps.
  3. Truncation semantics mismatch: in RK4, truncation(o) = trim(o, M) trims M
  TS representatives. In current TS Trotter, the same function trims M ordinary
  strings before TS rewrapping. That is a different and usually less favorable
  representation.

  In quick probes I saw the ordinary Trotter state hit the full cap while the
  TS-compressed representative count was much smaller, e.g. around 4096 ordinary
  terms versus roughly 1200–1400 TS representatives for the N=32, M=2^12
  example after a few steps. That is exactly the lost advantage.

  Why non-TS Trotter still uses the optimized path

  The non-TS method dispatch is src/evolution.jl:219-245.

  It calls:

  g = trotterize(H, dt; ...)
  trotter_step!(O, g; truncation=truncation)

  The optimized part is trotter_step! in src/trotter.jl:75-109.

  For each single-Pauli gate, it updates every Pauli string using only a cheap
  Pauli commutator:

  C, k = commutator(G, P)

  If G commutes with P, the string is unchanged. Otherwise the rotation gives a
  two-branch update:

  P -> cos(theta) P + sin(theta) P'

  This is much cheaper than RK4’s four full Hamiltonian commutators. In the
  non-TS case, both RK4 and Trotter work in the same ordinary Pauli-string
  representation, so Trotter’s lower constant wins.

  In the TS case, RK4 gets a TS-optimized algebra path, but Trotter falls back
  to ordinary full-gate evolution. That is the mismatch.

  Generic second-order TS Trotter design

  The generic fix is not to exploit MFIM-specific commutativity. Instead:

  ▎ Use the ordinary second-order Strang product over all translated Hamiltonian
  ▎ terms, but after each complete Trotter step, canonicalize the evolved 
  ▎ representative back into OperatorTS, apply TS-space truncation/dissipation, 
  ▎ and continue from the new representative.

  Mathematical structure

  Let

  H = \sum_{a=1}^m c_a \sum_{r \in G} T_r P_a T_r^\dagger

  where G is the translation group.

  The current method builds the ordinary list of all terms:

  h_j \in \{ c_a T_r P_a T_r^\dagger \}_{a,r}

  and uses Strang splitting:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2} L_{h_1}}
  e^{\frac{\Delta t}{2} L_{h_2}}
  \cdots
  e^{\Delta t L_{h_L}}
  \cdots
  e^{\frac{\Delta t}{2} L_{h_2}}
  e^{\frac{\Delta t}{2} L_{h_1}},

  where L_h(O) = i[h, O] / ħ in the Heisenberg picture.

  This is the standard second-order Strang/Trotter formula.

  The proposed TS method is:

  O_{n+1}^{TS}
  =
  \Pi_{TS}
  \left[
  S_2(\Delta t) \, \operatorname{rep}(O_n^{TS})
  \right],

  where Π_TS means canonical re-expression as an OperatorTS.

  This remains second order because:

  - S_2(dt) approximates the exact full flow with local error O(dt^3).
  - the exact flow under a translation-invariant Hamiltonian maps the TS
  subspace to itself.
  - canonical TS projection/re-expression after the complete symmetric step does
  not change the formal order on the TS subspace.

  Important: projection should happen after a whole Strang step, not after
  arbitrary individual gates, because single translated gates break translation
  symmetry.

  Proposed step algorithm

  Conceptually:

  function trotter_step_ts_projected(Hts, Ots, gates, dt; truncation,
  dissipation)
      # 1. Work from compressed representative.
      Or = representative(Ots)

      # 2. Apply a complete ordinary second-order Trotter step.
      trotter_step!(Or, gates; truncation=inner_truncation)

      # 3. Re-express result as translation-symmetric.
      Ots_next = OperatorTS{Ls,Ps}(Or)

      # 4. Apply dissipation/truncation in TS space, matching RK4 semantics.
      Ots_next = dissipation(Ots_next, dt)
      Ots_next = truncation(Ots_next)

      return Ots_next
  end

  Then the evolution loop stores Ots as the actual state, not Or.

  The next step begins again from:

  Or = representative(Ots)

  So the hot loop starts from the compressed representative every time.

  Gate construction

  The simplest version can still use:

  Hr = resum(H)
  gates = trotterize(Hr, dt; order=2)

  That is generic and already correct. But I would design the internal API so it
  does not conceptually require materializing Hr as the user-facing approach.

  Better internal representation:

  struct TrotterTSPlan
      translated_gates::Vector{TrotterGate}
      Ls
      Ps
  end

  or grouped as:

  struct TrotterOrbitGate
      generator::PauliStringTS
      coeff
  end

  with deterministic expansion into translated ordinary gates only when applying
  the complete Strang sweep.

  The first implementation can use the flat translated gates for simplicity; the
  design should leave room to avoid resum(H) allocation later.

  Expected performance

  Let:

  - m = length(Hts) = number of representative Hamiltonian terms.
  - Ntr = number of translations.
  - K = length(Ots) = number of TS representative strings after truncation.
  - Mfull = length(Or) = ordinary string count without TS recompression.

  Current TS Trotter roughly costs:

  O(2 * m * Ntr * Mfull)

  per second-order step.

  TS RK4 roughly costs:

  O(4 * m * Ntr * K)

  commutator work per step, plus intermediate allocations.

  Projected TS Trotter should cost closer to:

  O(2 * m * Ntr * K) + one TS canonicalization per step

  because every step starts from a TS-compressed representative.

  The constant per Trotter gate is also low: a single Pauli commutator and at
  most a two-branch update, from src/trotter.jl:75-109.

  So the expected generic win over RK4 is:

  1. about half as many Hamiltonian passes as RK4’s four stages,
  2. cheaper single-Pauli updates,
  3. TS-compressed state between steps,
  4. truncation in the same TS representation as RK4.

  This does not rely on commuting terms, 1D nearest-neighbor structure, or
  MFIM-specific layer decompositions.

  Optional later accelerations, not the main path

  Only after the generic projected Strang method is working and benchmarked, we
  can add optional structure-aware accelerations:

  - If translated copies of a Hamiltonian representative commute, apply the
  whole orbit as a commuting layer.
  - If terms admit graph coloring, apply color layers.
  - If the Hamiltonian has short-range locality, generate translated gates in
  cache-friendly geometric order.

  But these should be optimizations on top of the generic projected method, not
  the primary design.

  External references

  The product-formula basis is standard:

  - Trotter’s product formula: H. F. Trotter, “On the product of semi-groups of
  operators,” Proceedings of the AMS 10, 545–551 (1959). AMS link
  - Second-order symmetric splitting / Strang splitting: G. Strang, “On the
  Construction and Comparison of Difference Schemes,” SIAM Journal on Numerical 
  Analysis 5, 506–517 (1968). DOI
  - Higher-order Suzuki–Trotter background: N. Hatano and M. Suzuki, “Finding
  Exponential Product Formulas of Higher Orders,” in Quantum Annealing and Other
  Optimization Methods, Springer (2005). arXiv:math-ph/0506007
  - The current code’s single-Pauli rotation rule cites Pauli propagation ideas
  related to efficient Pauli-string rotation updates: Rudolph et al., “Pauli
  Propagation: A Computational Framework for Simulating Quantum Systems.”
  arXiv:2505.21606

  My recommendation: first implement/benchmark the
  project-after-each-complete-Strang-step design. It is generic, preserves
  second-order behavior, aligns truncation semantics with RK4, and directly
  targets the actual slowdown.

✻ Cooked for 6m 11s

※ recap: Goal: design a generic translation-symmetric second-order Trotter method that beats RK4 without special-case structure. Next action: discuss the proposed “project back to OperatorTS after each complete Strang step” design before implementing. (disable recaps in /config)

❯ So effectively, the TS operator is expanded, evolved throughout, and then     
  finally contracted? you are suggesting to apply the folding at each step      
  instead. Can we entirely bypass the expansion? That is, we apply the steps    
  directly in the efficient representation                                      

● Yes — roughly, but with one important correction:

  - H is literally expanded via resum(H).
  - O is not literally expanded as resum(O). The code evolves only
  representative(O).
  - However, after applying translated gates, that representative becomes a
  plain non-TS Operator, and it is kept that way throughout the step loop. It is
  only folded back into OperatorTS for saving/final return.

  So the practical effect is:

  TS H      -> expanded full H gates
  TS O      -> representative ordinary O
  evolution -> ordinary non-TS workspace
  output    -> folded back to OperatorTS only at save/final

  That is enough to lose much of the TS advantage.

  Can we bypass expansion entirely?

  Depends what “entirely” means.

  1. Bypass resum(H) materialization?

  Yes, definitely. This is generic.

  Instead of:

  Hr = resum(H)
  gates = trotterize(Hr, dt)

  we can build a lazy TS Trotter plan:

  for (c, p_ts) in H
      p = representative(p_ts)
      for shift in all_shifts(...)
          G = shift(p, Ls, Ps, shift)
          yield TrotterGate(G, theta)
      end
  end

  and then apply those generated gates without ever constructing the full Hr
  operator.

  This avoids the Hamiltonian expansion allocation and makes the codepath
  TS-aware, but it does not by itself remove the O(N) translated-gate work. A
  generic product formula over individual translated Pauli terms still has to
  visit those translated terms somehow.

  2. Keep O folded as OperatorTS between complete Trotter steps?

  Yes. This is the generic improvement I was suggesting.

  The safe generic algorithm is:

  Ots_n
    -> representative ordinary workspace
    -> apply one complete Strang sweep using lazy translated gates
    -> fold workspace back to OperatorTS
    -> apply truncation in TS space
    -> Ots_{n+1}

  So the state is unfolded only inside one complete Trotter step, then folded
  immediately. This is generic and should already fix the main slowdown because
  truncation is applied to TS representatives, not ordinary translated
  duplicates.

  3. Keep O as plain OperatorTS after every individual translated Pauli gate?

  Not generically correct.

  A single translated Pauli rotation, e.g.

  O \mapsto e^{i\theta G_r/2} O e^{-i\theta G_r/2},

  does not commute with translations unless G_r itself is translation invariant.
  It picks one origin/site/bond. Therefore it maps a TS operator to a generally
  non-TS operator.

  If we immediately fold after that one gate, we lose information about the
  relative position of the generated strings to the later gates in the ordered
  Trotter sweep.

  So this is only valid in special cases:

  - all translated copies commute,
  - the whole orbit is applied as an exactly translation-invariant primitive,
  - or we intentionally change the method to a projected/averaged product
  formula.

  Those may be useful, but they rely on additional structure or change the
  primitive.

  The best generic “direct efficient representation” design

  I think the right generic design is a TS-native composite step, not a
  TS-native single-gate step.

  Something like:

  struct TSTrotterPlan
      terms          # representative Hamiltonian terms and coefficients
      Ls
      Ps
      order
      hbar
  end

  Then:

  trotter_step!(Ots::OperatorTS, plan::TSTrotterPlan, dt; truncation)

  would internally:

  1. create a compact representative workspace:

  work = representative(Ots)
  2. generate translated gates lazily in Strang order;
  3. apply the existing optimized Pauli-rotation kernel to work;
  4. after the complete symmetric sweep, fold:

  Ots_next = OperatorTS{Ls,Ps}(work)
  5. apply truncation in TS space.

  This bypasses resum(H) and avoids carrying a fully expanded TS state across
  time steps. It still uses an ordinary/anchored workspace inside the step
  because generic individual gates break TS.

  That is probably the correct first implementation.

  Can we do even more directly than that?

  Possibly, but not with plain OperatorTS alone.

  To apply individual non-translation-invariant gates while still quotienting by
  translations, we would need an anchored TS workspace.

  Instead of storing just:

  PauliStringTS representative -> coeff

  we would store something like:

  (canonical PauliStringTS representative, frame/origin shift) -> coeff

  The extra frame/origin shift records where the current string sits relative to
  the ordered gate sweep. That lets us canonicalize strings without losing
  information needed by later translated gates.

  At the end of a complete translation-invariant/projection boundary, we can sum
  over frame shifts and return to ordinary OperatorTS.

  So the representation stack would be:

  OperatorTS               # efficient physical state
  AnchoredTSOperator        # internal non-TS-but-quotiented workspace
  OperatorTS               # after complete step

  This is more ambitious. It could reduce intermediate duplication further than
  a plain ordinary workspace, but it is more complex and should be justified by
  profiling after the simpler lazy/folded version.

  My refined recommendation

  Generic path, in order:

  1. Do not use resum(H) for Trotter planning. Generate translated gates lazily
  from H::OperatorTS.
  2. Keep the public/state-level evolution in OperatorTS. Fold after every
  complete second-order Trotter step.
  3. Apply truncation after folding, so trim(o, M) means “keep M TS
  representatives,” matching RK4 semantics.
  4. Only if profiling still shows within-step ordinary workspace growth
  dominates, design an AnchoredTSOperator internal workspace.
  5. Treat commuting-orbit layers, even/odd coloring, etc. as optional later
  accelerations, not the core method.

  So: yes, we can bypass Hamiltonian expansion entirely and make the Trotter 
  step TS-native at the step boundary. But no, a generic method cannot apply 
  each individual translated Pauli rotation purely inside OperatorTS without 
  extra anchored information, because individual gates break translation
  symmetry.

✻ Churned for 2m 37s

❯ I see. Can we apply the entire trotter step in one go? That is, we do not     
  expand resum (even lazily), and instead apply the TS preserving trotter step  
  directly in the TS subspace?                                                  
● Short answer: not generically with only the current OperatorTS representation 
  and the single-Pauli Trotter primitive.

  We can define a TS-preserving “whole Trotter step” mathematically, but
  applying it efficiently and generically requires either:

  1. an internal non-TS/anchored workspace,
  2. a reduced TS superoperator/Krylov-style method,
  3. or extra structure such as commuting translated copies.

  So the practical generic path is probably:

  avoid materialized resum(H)
  generate translated gates internally/lazily or through an anchored quotient
  apply one complete Strang sweep
  fold once per complete step
  truncate in TS space

  But a true “apply directly inside plain OperatorTS at all times” is not
  generally available.

  ---
  Why plain OperatorTS is not enough

  OperatorTS stores translation orbits:

  representative Pauli string -> coefficient

  It intentionally forgets the absolute origin.

  That is fine for operations that are themselves translation-equivariant, like:

  commutator(Hts, Ots)

  because the TS commutator can be written as a sum over relative shifts and
  then folded back into representatives.

  But a single translated Trotter gate,

  O \mapsto e^{i \theta G_r/2} O e^{-i \theta G_r/2},

  is not translation-equivariant. It singles out site/bond/origin r.

  After applying such a gate, the intermediate operator generally has an
  origin-dependent defect. If we immediately fold it into OperatorTS, we lose
  information needed by the next translated gate in the ordered Trotter sweep.

  A toy generic example is an orbit term like:

  H_a = \sum_i X_i Z_{i+1}.

  Translated copies of X_i Z_{i+1} need not commute with neighboring copies. A
  sweep over i = 1, 2, ... creates intermediate strings whose relation to the
  current sweep position matters. Plain OperatorTS cannot encode “this string is
  three sites behind the current gate” because it quotients out the origin.

  So the problem is not just implementation detail; it is a representation
  issue.

  ---
  What about applying the entire step, not prefixes?

  There are two possible meanings.

  Meaning 1: Product over all translated Pauli gates, then fold

  This is closest to the current Trotter method.

  For a full expanded Hamiltonian:

  H = \sum_{a,r} h_{a,r},

  ordinary second-order Strang uses something like:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2}L_{h_1}}
  e^{\frac{\Delta t}{2}L_{h_2}}
  \cdots
  e^{\Delta t L_{h_L}}
  \cdots
  e^{\frac{\Delta t}{2}L_{h_2}}
  e^{\frac{\Delta t}{2}L_{h_1}}.

  This full composite map can be applied to a representative and folded after
  the complete sweep. That is the “fold every step” design.

  However, evaluating the composite map still requires tracking the non-TS
  intermediate state somehow. If we do not use ordinary Operator, we need a
  richer internal representation, e.g. an anchored TS workspace.

  So this can avoid resum(H) materialization, but it cannot generally avoid the
  logical translated-gate sweep unless we change the numerical method.

  Meaning 2: Split by TS orbit Hamiltonians

  We could define a different, very natural TS-preserving Strang formula:

  H = \sum_{a=1}^m H_a^{TS},
  \qquad
  H_a^{TS} = c_a \sum_r T_r P_a T_r^\dagger.

  Then use:

  S_2^{TS}(\Delta t)
  =
  e^{\frac{\Delta t}{2}L_{H_1^{TS}}}
  e^{\frac{\Delta t}{2}L_{H_2^{TS}}}
  \cdots
  e^{\Delta t L_{H_m^{TS}}}
  \cdots
  e^{\frac{\Delta t}{2}L_{H_2^{TS}}}
  e^{\frac{\Delta t}{2}L_{H_1^{TS}}}.

  This is genuinely TS-preserving at every factor. It also does not expand
  resum(H) in the formula.

  But the hard part is applying one factor:

  e^{\tau L_{H_a^{TS}}} O
  =
  e^{i \tau H_a^{TS}} O e^{-i \tau H_a^{TS}}.

  H_a^{TS} is not a single Pauli string. It is a sum over translated Pauli
  strings. If those translated strings do not all commute, then there is no
  cheap two-branch Pauli-rotation update.

  The existing fast Trotter kernel works because for a single Pauli generator G:

  e^{i\theta G/2} P e^{-i\theta G/2}

  is either unchanged or becomes a two-term combination. For a generic orbit sum
  H_a^{TS}, the exponential can generate many nested commutators. Applying it
  exactly is basically as hard as solving a smaller but still many-body
  evolution problem.

  So this second meaning is mathematically elegant but not automatically faster.
  It turns the problem into “how do we exponentiate an orbit Hamiltonian in TS
  space?” Generic answers are RK, Krylov/Lanczos, Chebyshev, or an inner Trotter
  over translations — which brings us back to translated gates.

  ---
  Is the current expanded Strang step itself TS-preserving?

  Important subtlety: an ordered product over translated local gates is not
  necessarily exactly translation-invariant at finite dt.

  If:

  U = e^{A_1} e^{A_2} \cdots e^{A_N},

  then translating gives a cyclically/permuted product:

  T U T^\dagger = e^{A_2} e^{A_3} \cdots e^{A_1}.

  If the A_i do not commute, this is generally not equal to U.

  A palindromic Strang order gives second-order time symmetry, but it does not
  automatically restore exact translation symmetry of the finite-step map.

  The exact flow under H is TS-preserving. The Trotterized finite-step map may
  only preserve TS up to its product-formula error unless the ordering/layers
  are chosen with symmetry in mind or the relevant terms commute.

  This is another reason a “direct TS-subspace product step” is not automatic:
  we first need to define which TS-preserving approximation we want.

  ---
  Generic options

  Option A: Lazy translated gates + fold after each complete step

  This is the pragmatic generic fix.

  Do not build:

  Hr = resum(H)

  Instead build a TSTrotterPlan from representative TS terms and translations.

  But internally, for each complete Strang step, still apply the translated
  single-Pauli rotations in the right order, then fold:

  Ots_n
    -> representative workspace
    -> complete Strang sweep over translated gates
    -> OperatorTS fold
    -> truncation in TS space
    -> Ots_{n+1}

  Pros:

  - generic;
  - preserves existing Trotter semantics;
  - avoids materialized resum(H);
  - restores TS compression between steps;
  - matches RK4 truncation semantics better.

  Cons:

  - still has O(number of translations) gate work per step;
  - still needs a non-TS intermediate workspace inside the step.

  This is what I would implement first.

  ---
  Option B: Anchored TS workspace

  To get closer to “direct efficient representation,” introduce an internal
  representation like:

  (canonical TS representative, anchor/frame/origin) -> coefficient

  This keeps enough origin information to apply an ordered translated-gate sweep
  without expanding all strings as ordinary PauliStrings.

  Conceptually:

  OperatorTS
    -> AnchoredTSOperator
    -> apply translated gates while updating anchors
    -> fold anchors away
    -> OperatorTS

  Pros:

  - generic;
  - avoids plain full ordinary expansion of the state;
  - preserves information that plain OperatorTS loses.

  Cons:

  - more design work;
  - more complicated hashing/canonicalization;
  - may or may not beat the simpler fold-per-step version enough to justify
  itself.

  This is the generic route if profiling shows within-step ordinary workspace
  growth still dominates.

  ---
  Option C: Reduced TS superoperator

  For fixed H, dt, and a fixed truncated TS basis, one could precompute the full
  step as a matrix/sparse map:

  S_{TS}: \mathcal{B}_{TS} \to \mathcal{B}_{TS}.

  Then each time step is just a TS-space matrix-vector multiply.

  Pros:

  - applies the whole step “in one go” after precomputation;
  - directly TS-subspace based.

  Cons:

  - basis can be huge;
  - trim(o, M) is adaptive, so the active basis changes;
  - building the superoperator columns still requires applying the step somehow;
  - memory likely bad for large benchmarks.

  Potentially useful for many repeated runs on the same fixed basis, but not my
  first generic recommendation.

  ---
  Option D: Orbit-level TS Strang with exact orbit exponentials

  This is the cleanest mathematically:

  e^{\tau L_{H_a^{TS}}}

  directly in TS space.

  But generic exact orbit exponentials are not cheap unless the translated terms
  commute or have special solvable structure. Otherwise this becomes
  RK/Krylov-like evolution under H_a^{TS}.

  So this is not a generic Trotter-speed solution.

  ---
  Conclusion

  If the question is:

  ▎ Can we bypass materializing resum(H)?

  Yes.

  If the question is:

  ▎ Can we avoid even lazy translated-gate application and apply the whole 
  ▎ generic Trotter step directly in plain OperatorTS?

  Not generally, not while retaining the fast single-Pauli Trotter mechanism.
  Plain OperatorTS has intentionally discarded origin information, and generic
  ordered Trotter prefixes need that information.

  The most defensible generic design is therefore:

  TSTrotterPlan:
      stores representative TS Hamiltonian terms
      generates translated gates internally, not via resum(H)

  trotter_step_ts!:
      starts from OperatorTS
      uses a representative or anchored workspace for one complete Strang step
      folds back to OperatorTS after the complete step
      truncates in TS space

  Then, if that still is not enough, the next generic upgrade is an anchored TS 
  workspace, not a plain OperatorTS-only step.

  References: the product formula foundation is Trotter’s original semigroup
  product formula, Strang’s second-order symmetric splitting, and
  Suzuki/Hatano’s higher-order product-formula framework: Trotter 1959, Strang 
  1968, Hatano & Suzuki 2005.

✻ Sautéed for 2m 34s

❯ Lets go with the first option. Implement the lazy generation of resum         
  elements                                                                      

  Searched for 2 patterns, read 3 files (ctrl+o to expand)

● Update(src/trotter.jl)
  ⎿  Added 64 lines, removed 1 line
       42      return gates
       43  end
       44  
       45 +_ordinary_paulistringtype(::Type{<:PauliStringTS{Ls,Ps,U}}) wh
          +ere {Ls,Ps,U} = PauliString{Base.prod(Ls),U}                  
       46 +                                                              
       47  """
       48 +    TrotterGatesTS                                            
       49 +                                                              
       50 +Lazy first- or second-order Trotter gate list for a translatio
          +n-symmetric Hamiltonian.                                      
       51 +The generators are the translated representatives of `H`; they
          + are produced on demand                                       
       52 +instead of materializing `resum(H)`.                          
       53 +"""                                                           
       54 +struct TrotterGatesTS{P<:PauliString,O<:Operator{<:PauliString
          +TS},S,D<:Real,R<:Real} <: AbstractVector{TrotterGate{P,Float64
          +}}                                                            
       55 +    H::O                                                      
       56 +    shifts::Vector{S}                                         
       57 +    dt::D                                                     
       58 +    hbar::R                                                   
       59 +    order::Int                                                
       60 +    heisenberg::Bool                                          
       61 +end                                                           
       62 +                                                              
       63 +Base.IndexStyle(::Type{<:TrotterGatesTS}) = IndexLinear()     
       64 +Base.size(g::TrotterGatesTS) = (length(g),)                   
       65 +                                                              
       66 +function Base.length(g::TrotterGatesTS)                       
       67 +    L = length(g.H) * length(g.shifts)                        
       68 +    return (g.order == 1 || L <= 1) ? L : 2L - 1              
       69 +end                                                           
       70 +                                                              
       71 +function _expanded_trotter_gate(g::TrotterGatesTS{P}, i::Int, 
          +scale::Real) where {P}                                        
       72 +    L = length(g.H)                                           
       73 +    shift_index = (i - 1) ÷ L + 1                             
       74 +    term_index = (i - 1) % L + 1                              
       75 +    c, p = g.H[term_index]                                    
       76 +    generator = shift(representative(p), qubitsize(g.H), perio
          +dicflags(g.H), g.shifts[shift_index])                         
       77 +    theta = Float64(scale * _trotter_theta(c, g.dt, g.hbar, g.
          +heisenberg))                                                  
       78 +    return TrotterGate{P,Float64}(generator, theta)           
       79 +end                                                           
       80 +                                                              
       81 +function Base.getindex(g::TrotterGatesTS, i::Int)             
       82 +    @boundscheck checkbounds(g, i)                            
       83 +    L = length(g.H) * length(g.shifts)                        
       84 +    if g.order == 1 || L <= 1                                 
       85 +        return _expanded_trotter_gate(g, i, 1)                
       86 +    elseif i < L                                              
       87 +        return _expanded_trotter_gate(g, i, 1 / 2)            
       88 +    elseif i == L                                             
       89 +        return _expanded_trotter_gate(g, L, 1)                
       90 +    else                                                      
       91 +        return _expanded_trotter_gate(g, 2L - i, 1 / 2)       
       92 +    end                                                       
       93 +end                                                           
       94 +                                                              
       95 +function _trotter_gates_ts(H::Operator{<:PauliStringTS}, dt::R
          +eal, hbar::Real, heisenberg::Bool, order::Integer)            
       96 +    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 
          +2, got $order"))                                              
       97 +    norm(H - H') > 1e-10 && throw(ArgumentError("Hamiltonian m
          +ust be Hermitian for Trotter splitting"))                     
       98 +    P = _ordinary_paulistringtype(paulistringtype(H))         
       99 +    shifts = collect(all_shifts(paulistringtype(H)))          
      100 +    return TrotterGatesTS{P,typeof(H),eltype(shifts),typeof(dt
          +),typeof(hbar)}(H, shifts, dt, hbar, Int(order), heisenberg)  
      101 +end                                                           
      102 +                                                              
      103 +"""                                                           
      104      trotterize(H::Operator, dt::Real; order=2, heisenberg=true
           , hbar=1)
      105  
      106  Build a first-order (`order=1`, Lie) or second-order (`order=2
           `, Strang) Trotter list that approximates
      107  `exp(im * H * dt / hbar)` (Heisenberg) or the conjugate sequen
           ce (Schrödinger / density matrix).
      108  Each gate uses [`pauli_rotation`](@ref) with the returned `the
           ta` field.
      109  
       52 -For `H::Operator{<:PauliStringTS}`, see the specialized [`trot
          -terize`](@ref) that calls [`resum`](@ref) first.              
      110 +For `H::Operator{<:PauliStringTS}`, the returned gate list laz
          +ily generates translated                                      
      111 +representative terms instead of materializing `resum(H)`.     
      112  """
      113  function trotterize(H::Operator, dt::Real; order::Integer=2, h
           eisenberg::Bool=true, hbar::Real=1)
      114      order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 
           2, got $order"))
     ...
      124      end
      125  end
      126  
      127 +function trotterize(H::Operator{<:PauliStringTS}, dt::Real; or
          +der::Integer=2, heisenberg::Bool=true, hbar::Real=1)          
      128 +    return _trotter_gates_ts(H, dt, hbar, heisenberg, order)  
      129 +end                                                           
      130 +                                                              
      131  """
      132  trotter_step!(O::Operator, gates::AbstractVector{<:TrotterGate
           }; truncation::Function=identity, truncate_every::Int=1)
      133  

● Update(src/evolution.jl)
  ⎿  Added 15 lines, removed 12 lines
      252      Ps = periodicflags(O)
      253      n = length(tspan)
      254      history = _alloc_history(fout, O, n)
      255 -    # Evolve the representative as a dense `Operator`: the gat
          -es are built from the                                         
      256 -    # resummed (translation-unfolded) `H`, so applying them to
          - the representative term                                      
      257 -    # reproduces the translation-symmetric dynamics. `represen
          -tative` returns a fresh                                       
      258 -    # `Operator`, so `trotter_step!`'s in-place mutation never
          - aliases the caller's `O`.                                    
      259 -    Or = representative(O)                                    
      260 -    Hr = resum(H)                                             
      255 +    # Keep the state in the translation-symmetric representati
          +on between full                                               
      256 +    # Trotter steps. Individual translated Pauli rotations are
          + not themselves                                               
      257 +    # translation symmetric, so each complete step uses a repr
          +esentative workspace                                          
      258 +    # and folds it back afterwards.                           
      259 +    Ots = copy(O)                                             
      260  
      261      # Cache gates when the save spacing is uniform; rebuild pe
           r step otherwise.
      262 +    # `trotterize(H, dt)` for `H::OperatorTS` generates transl
          +ated gates lazily;                                            
      263 +    # it does not materialize `resum(H)`.                     
      264      dt0 = n > 1 ? (tspan[2] - tspan[1]) : zero(eltype(tspan))
      265      uniform = n > 1 && all(i -> tspan[i + 1] - tspan[i] ≈ dt0,
            1:(n - 1))
      266      gates_cached = method.gates !== nothing ? method.gates :
      267                     (uniform && n > 1 ?
      267 -                    trotterize(Hr, dt0; order=method.order, he
          -isenberg=true, hbar=hbar) :                                   
      268 +                    trotterize(H, dt0; order=method.order, hei
          +senberg=true, hbar=hbar) :                                    
      269                      nothing)
      270  
      271      for i in ProgressBar(1:(n - 1))
      272          dt = tspan[i + 1] - tspan[i]
      273          g = gates_cached !== nothing ? gates_cached :
      273 -            trotterize(Hr, dt; order=method.order, heisenberg=
          -true, hbar=hbar)                                              
      274 +            trotterize(H, dt; order=method.order, heisenberg=t
          +rue, hbar=hbar)                                               
      275 +        Or = representative(Ots)                              
      276          trotter_step!(Or, g; truncation=truncation)
      275 -        Or = dissipation(Or, dt)                              
      276 -        Or = truncation(Or)                                   
      277 -        _save!(history, fout, OperatorTS{Ls,Ps}(Or), i + 1)   
      277 +        Ots = OperatorTS{Ls,Ps}(Or)                           
      278 +        Ots = dissipation(Ots, dt)                            
      279 +        Ots = truncation(Ots)                                 
      280 +        _save!(history, fout, Ots, i + 1)                     
      281      end
      279 -    return EvolutionResult(collect(tspan), history, OperatorTS
          -{Ls,Ps}(Or))                                                  
      282 +    return EvolutionResult(collect(tspan), history, Ots)      
      283  end
      284  
      285  function _evolve(::Trotter, H::AbstractOperator, O::AbstractOp
           erator, tspan;

  Read 2 files, ran 2 shell commands (ctrl+o to expand)

● Update(src/trotter.jl)
  ⎿  Added 2 lines, removed 1 line
       97      norm(H - H') > 1e-10 && throw(ArgumentError("Hamiltonian m
           ust be Hermitian for Trotter splitting"))
       98      P = _ordinary_paulistringtype(paulistringtype(H))
       99      shifts = collect(all_shifts(paulistringtype(H)))
      100 -    return TrotterGatesTS{P,typeof(H),eltype(shifts),typeof(dt
          -),typeof(hbar)}(H, shifts, dt, hbar, Int(order), heisenberg)  
      100 +    Hcopy = copy(H)                                           
      101 +    return TrotterGatesTS{P,typeof(Hcopy),eltype(shifts),typeo
          +f(dt),typeof(hbar)}(Hcopy, shifts, dt, hbar, Int(order), heise
          +nberg)                                                        
      102  end
      103  
      104  """

  Read 1 file, ran 1 shell command (ctrl+o to expand)

● Update(test/trotter.jl)
  ⎿  Added 6 lines
      40      gh = trotterize(H1, 0.11; heisenberg=true)
      41      gs = trotterize(H1, 0.11; heisenberg=false)
      42      @test gh[1].theta == -gs[1].theta
      43 +                                                               
      44 +    Hts = OperatorTS{(4,)}(H)                                  
      45 +    gts = trotterize(Hts, 0.1; order=2, heisenberg=true)       
      46 +    @test gts isa AbstractVector{<:TrotterGate}                
      47 +    @test !(gts isa Vector)                                    
      48 +    @test length(gts) == length(trotterize(resum(Hts), 0.1; ord
         +er=2, heisenberg=true))                                        
      49  end
      50  
      51  @testset "trotter_step! - commuting terms vs exact Heisenberg" 
          begin

● Update(test/evolution.jl)
  ⎿  Added 2 lines
      32      res3 = evolve(H, O0, times; method = DOPRI5(), fout=fout, t
          runcation=truncation).history
      33      res4 = evolve(resum(H), resum(O0), times; method = Trotter(
          ), fout=fout_trotter, truncation=truncation).history
      34      res5 = evolve(H, O0, times; method = Trotter(), fout=fout, 
          truncation=truncation).history
      35 +    res6 = evolve(H, O0, [0, dt]; method = Trotter(), fout=fout
         +, truncation=truncation).final                                 
      36  
      37      @test norm(res2 .- res1)/norm(res1) < 1e-7
      38      @test norm(res3 .- res1)/norm(res1) < 1e-7
      39      @test norm(res4 .- res1)/norm(res1) < 1e-5
      40      @test norm(res5 .- res1)/norm(res1) < 1e-4
      41 +    @test res6 isa typeof(O0)                                  
      42  
      43      dissipation(O, dt) = add_noise(O, 0.1*dt)
      44      res1 = evolve(H, O0, times; method = RK4(), fout=fout, diss
          ipation=dissipation).history

  Ran 1 shell command (ctrl+o to expand)

● Update(test/trotter.jl)
  ⎿  Added 4 lines, removed 1 line
      41      gs = trotterize(H1, 0.11; heisenberg=false)
      42      @test gh[1].theta == -gs[1].theta
      43  
      44 -    Hts = OperatorTS{(4,)}(H)                                  
      44 +    Hts0 = Operator(4)                                         
      45 +    Hts0 += 0.5, "X", 1                                        
      46 +    Hts0 += 0.3, "Z", 2                                        
      47 +    Hts = OperatorTS{(4,)}(Hts0)                               
      48      gts = trotterize(Hts, 0.1; order=2, heisenberg=true)
      49      @test gts isa AbstractVector{<:TrotterGate}
      50      @test !(gts isa Vector)

  Read 1 file, ran 2 shell commands (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 2 lines
      16  when `tspan` is uniformly spaced.
      17  
      18  Implemented for `H::Operator`/`O::Operator` and for `H::Operato
          rTS`/`O::OperatorTS`. In
      19 -the translation-symmetric case the gates are built from `resum(
         -H)` and applied to the                                         
      20 -representative term, so a precomputed `gates` list must come fr
         -om `trotterize(resum(H), dt)`.                                 
      19 +the translation-symmetric case the gate list lazily generates t
         +he translated representatives                                  
      20 +of `H` instead of materializing `resum(H)`, applies each comple
         +te step to a representative                                    
      21 +workspace, then folds the result back to `OperatorTS` before tr
         +uncation and output.                                           
      22  """
      23  struct Trotter <: AbstractEvolutionMethod
      24      order::Int

  Ran 1 shell command (ctrl+o to expand)
● Implemented the first option: lazy generation of resum(H) elements for TS
  Trotter.

  What changed

  - Added a lazy TrotterGatesTS type in src/trotter.jl.
    - It subtypes AbstractVector{TrotterGate{P,Float64}}, so existing
  trotter_step! can consume it unchanged.
    - It stores:
        - copied TS Hamiltonian H,
      - translation shifts,
      - dt, hbar, order, heisenberg.
    - getindex generates each translated ordinary gate on demand via:

  shift(representative(p), qubitsize(H), periodicflags(H), shift)
    - It no longer materializes Hr = resum(H) just to build gates.
  - Added specialized:

  trotterize(H::Operator{<:PauliStringTS}, dt; ...)

  - which returns TrotterGatesTS.
  - Updated TS evolve(...; method=Trotter()) in src/evolution.jl.
    - Removed:

  Hr = resum(H)
  trotterize(Hr, ...)
    - Replaced with:

  trotterize(H, ...)
    - Keeps O as OperatorTS between complete Trotter steps:
        i. Or = representative(Ots)
      ii. apply full Trotter step to Or
      iii. fold back with OperatorTS{Ls,Ps}(Or)
      iv. apply dissipation and truncation in TS space
  - Updated docs in src/evolution.jl to reflect this new behavior.
  - Added tests:
    - trotterize(Hts, ...) returns a lazy non-Vector AbstractVector of gates.
    - Lazy TS gate list has the same length as trotterize(resum(Hts), ...).
    - TS Trotter final state remains the original OperatorTS type.

  Validation

  Passed:

  test/trotter.jl
  test/evolution.jl

  via:

  julia --project=/home/dhruva/Programs/unitaryHack/PauliStrings.jl -e 'using 
  PauliStrings; include("test/trotter.jl"); include("test/evolution.jl")'

  Performance note

  I also ran a quick benchmark after implementation. Results were not yet faster
  than RK4:

  N=32, M=12, 50 steps: RK4 ≈ 1.32s, Trotter ≈ 1.81s
  N=64, M=14, 20 steps: RK4 ≈ 2.35s, Trotter ≈ 4.47s

  So this successfully removes materialized resum(H) and fixes the state-level
  TS folding semantics, but it exposes the next bottleneck:
  folding/canonicalizing every step plus still applying all translated gates
  individually. The next generic optimization would likely be improving the
  internal workspace/folding cost, e.g. an anchored TS workspace or more
  targeted folding, rather than relying on MFIM-specific commuting structure.

✻ Brewed for 10m 3s

❯ Can you check the time for the same tests before and after the lazy resum?    
  You can git stash, run, and git stash apply                                   

● Write(~/.claude/jobs/513a65b3/tmp/bench_ts_trotter.jl)
  ⎿  Wrote 62 lines to
     ../../../.claude/jobs/513a65b3/tmp/bench_ts_trotter.jl
      1 using PauliStrings
      2 
      3 function MFIM(N, h)
      4     H = Operator(N)
      5     H += -h, "X", 1
      6     H += -h/2, "Z", 2
      7     H += "Z", 1, "Z", 2
      8     return OperatorTS{(N,)}(H)
      9 end
     … +59 lines (ctrl+o to expand)

  Ran 3 shell commands (ctrl+o to expand)

● I checked before/after using the same benchmark script and git stash as
  requested.

  Benchmark script used:

  - MFIM TS Hamiltonian
  - dt = 0.1
  - truncation(o) = trim(o, 2^M)
  - repeated 3 times per case after warmup
  - reported best time

  Results

  ┌────────────────┬────────────────────┬─────────┬───────────┬─────────────┐
  │      Case      │      Version       │  RK4    │ Trotter   │ Trotter/RK4 │
  │                │                    │  best   │   best    │             │
  ├────────────────┼────────────────────┼─────────┼───────────┼─────────────┤
  │ N=32, M=12,    │  before lazy resum │  1.173s │    1.342s │       1.14x │
  │ steps=50       │                    │         │           │             │
  ├────────────────┼────────────────────┼─────────┼───────────┼─────────────┤
  │ N=32, M=12,    │ after lazy resum + │  1.169s │    1.690s │       1.44x │
  │ steps=50       │          fold/step │         │           │             │
  ├────────────────┼────────────────────┼─────────┼───────────┼─────────────┤
  │ N=64, M=14,    │  before lazy resum │  2.285s │    4.078s │       1.78x │
  │ steps=20       │                    │         │           │             │
  ├────────────────┼────────────────────┼─────────┼───────────┼─────────────┤
  │ N=64, M=14,    │ after lazy resum + │  2.231s │    4.254s │       1.91x │
  │ steps=20       │          fold/step │         │           │             │
  └────────────────┴────────────────────┴─────────┴───────────┴─────────────┘

  Interpretation

  The current implementation I added is slower than baseline for these tests.

  Why:

  - Avoiding materialized resum(H) saves only the one-time gate-list
  construction cost.
  - The dominant cost is still applying all translated gates every step.
  - The new version additionally folds back to OperatorTS after every complete
  Trotter step and truncates in TS space, which adds canonicalization overhead.
  - For these sizes, that added per-step folding cost is larger than the saved
  resum(H) materialization.

  One caveat: the “before” script reports gates=5 because trotterize(H, dt) was
  not specialized in baseline and returns only representative TS gates when
  called directly. But the actual baseline evolve(H, O0; method=Trotter()) path
  still used Hr = resum(H) internally, as seen in src/evolution.jl, so the
  timing comparison is valid.

  So: lazy resum generation alone is not the bottleneck fix. The next thing to
  target should be reducing either:

  1. the repeated folding/canonicalization cost, or
  2. the per-step translated-gate sweep cost.

  The working tree has been restored with the lazy changes reapplied.

✻ Crunched for 2m 14s

❯ your primary hypothesis was that the major slow down was happening due to the 
  expansion of the TS variant to the normal variant. If folding is not          
  providing a major spedup then the hypothesis is wrong. Instead of moving      
  forward trying to only reduce the run time, actually figure out where the     
  slow down is.                                                                 

● You're right. The timing result falsifies the “materialized resum(H)/state
  folding is the dominant cost” hypothesis for this benchmark. I’ll stop
  optimizing that path and profile where the time is actually going.

● Write(~/.claude/jobs/513a65b3/tmp/profile_ts_components.jl)
  ⎿  Wrote 145 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_ts_components.jl
      1 using PauliStrings
      2 
      3 function MFIM(N, h)
      4     H = Operator(N)
      5     H += -h, "X", 1
      6     H += -h/2, "Z", 2
      7     H += "Z", 1, "Z", 2
      8     return OperatorTS{(N,)}(H)
      9 end
     … +150 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● The component profile already points away from resum(H) materialization, but
  I’m going one level deeper before concluding. The first clear signal is that
  Trotter is applying/truncating per gate — hundreds of times per step — while
  RK4 truncates only a handful of times per step and uses the TS commutator
  kernel. I’ll quantify that now.

● Write(~/.claude/jobs/513a65b3/tmp/profile_trotter_truncate_every.jl)
  ⎿  Wrote 93 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_trotter_truncate_every.jl
      1 using PauliStrings
      2 
      3 function MFIM(N, h)
      4     H = Operator(N)
      5     H += -h, "X", 1
      6     H += -h/2, "Z", 2
      7     H += "Z", 1, "Z", 2
      8     return OperatorTS{(N,)}(H)
      9 end
     … +95 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

●

  Read 1 file, ran 2 shell commands (ctrl+o to expand)

● Background command "Profile Trotter truncation frequency impact" was stopped

● Write(~/.claude/jobs/513a65b3/tmp/profile_rk4_internals.jl)
  ⎿  Wrote 48 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_rk4_internals.jl
       1 using PauliStrings
       2 
       3 function MFIM(N, h)
       4     H = Operator(N)
       5     H += -h, "X", 1
       6     H += -h/2, "Z", 2
       7     H += "Z", 1, "Z", 2
       8     return OperatorTS{(N,)}(H)
       9 end
      10 function Xtot(N)
     … +43 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)
● You’re right: the “expanded TS state is the major slowdown” hypothesis was
  wrong, at least for this benchmark. I profiled components instead of
  continuing down that path.

  Corrected finding

  The dominant slowdown is not one-time resum(H) materialization, and folding
  every step does not help. The actual bottleneck is:

  ▎ Current Trotter does hundreds of sequential single-Pauli gate updates per 
  ▎ time step, and it truncates after every gate by default. RK4 does only four 
  ▎ TS commutators and ~five truncations per step.

  For MFIM:

  length(Hts) = 3
  length(resum(Hts)) = 3N
  Strang gates = 2 * 3N - 1 = 6N - 1

  So:

  - N=32: 191 Trotter gates per time step
  - N=64: 383 Trotter gates per time step

  And trotter_step! currently does:

  for each gate
      update all strings
      truncation(O2)   # default truncate_every = 1
  end

  So at N=64, each physical time step causes 383 calls to trim, not one.

  RK4, by contrast, does four TS commutator evaluations and calls truncation
  around five times per step.

  Component profile

  I profiled short runs to avoid misleading long-runtime effects.

  N=32, M=12, steps=10

  RK4

  total       0.259s
  commutator  0.103s  over 40 calls
  trim        0.139s  over 50 calls

  Baseline current Trotter: resum(H) once, ordinary workspace across steps

  total       0.308s
  resum(H)    0.054s  one-time
  gates build 0.113s  one-time
  trim        0.087s  over 1920 calls
  gates       191

  Important: the actual loop time is 0.308s; one-time resum + gate construction
  is not the main loop cost.

  Lazy/fold-every-step Trotter

  total       0.227s
  lazy gates  0.005s  one-time
  step body   0.224s
  fold        0.002s
  trim        0.078s over 1920 calls

  For this short case lazy was faster than baseline, but this did not hold
  robustly for longer/larger cases because folding changes state/truncation
  dynamics and can increase work later.

  N=64, M=14, steps=8

  RK4

  total       0.539s
  commutator  0.453s over 32 calls
  trim        0.011s over 40 calls

  Baseline current Trotter

  total       1.503s
  resum(H)    0.053s one-time
  gates build 0.110s one-time
  trim        0.368s over 3072 calls
  gates       383

  Lazy/fold-every-step Trotter

  total       1.909s
  lazy gates  0.006s one-time
  step body   1.892s
  fold        0.017s
  trim        0.587s over 3072 calls

  Here the cost is clearly not resum(H). It is the gate sweep plus the per-gate
  trimming.

  Truncation frequency experiment

  I then varied truncate_every in trotter_step!.

  N=32, M=12, steps=5, baseline Trotter

  ┌────────────────┬────────┬────────────┬─────────────────────────┐
  │ truncate_every │ total  │ trim calls │ max intermediate length │
  ├────────────────┼────────┼────────────┼─────────────────────────┤
  │              1 │ 0.284s │        960 │                   5,289 │
  ├────────────────┼────────┼────────────┼─────────────────────────┤
  │              4 │ 0.094s │        240 │                   5,893 │
  ├────────────────┼────────┼────────────┼─────────────────────────┤
  │             16 │ 0.088s │         60 │                  10,934 │
  ├────────────────┼────────┼────────────┼─────────────────────────┤
  │             64 │ 0.509s │         15 │                 182,329 │
  ├────────────────┼────────┼────────────┼─────────────────────────┤
  │            191 │ 12.58s │         10 │               5,918,621 │
  └────────────────┴────────┴────────────┴─────────────────────────┘

  N=64, M=14, steps=4, baseline Trotter

  ┌────────────────┬─────────┬────────────┬─────────────────────────┐
  │ truncate_every │  total  │ trim calls │ max intermediate length │
  ├────────────────┼─────────┼────────────┼─────────────────────────┤
  │              1 │  0.699s │       1536 │                  23,683 │
  ├────────────────┼─────────┼────────────┼─────────────────────────┤
  │              4 │  0.441s │        384 │                  23,730 │
  ├────────────────┼─────────┼────────────┼─────────────────────────┤
  │             16 │  0.422s │         96 │                  35,907 │
  ├────────────────┼─────────┼────────────┼─────────────────────────┤
  │             64 │  0.704s │         24 │                  67,257 │
  ├────────────────┼─────────┼────────────┼─────────────────────────┤
  │            383 │ 173.90s │          8 │              19,715,254 │
  └────────────────┴─────────┴────────────┴─────────────────────────┘

  This is very informative:

  - Truncating every gate is too frequent.
  - Truncating only at the end is catastrophically bad because the branching
  explosion is enormous.
  - There is a sweet spot around truncate_every = 4–16 for these cases.

  So Trotter is losing because it is caught between two bad regimes:

  1. trim too frequently → many expensive sort/partialsort calls,
  2. trim too rarely → intermediate operator explodes.

  The current default truncate_every=1 chooses safety but pays a huge overhead.

  Why RK4 wins in TS

  RK4’s TS commutator path is actually quite good.

  For N=64, M=14, steps=8, RK4 spent:

  0.448s in commutators
  0.033s in trim
  0.051s in operator additions
  0.532s total

  The TS binary kernel does one relative-translation loop per Hamiltonian
  representative pair and folds into PauliStringTS immediately. It performs few
  high-level steps and few truncations.

  Trotter performs many more high-level sequential updates:

  383 gate updates per physical step
  383 truncation opportunities per physical step

  Even though each gate update is cheap, there are so many gates that the total
  cost is worse.

  Revised diagnosis

  The slowdown is caused by the granularity mismatch:

  - Non-TS Trotter is good because the number of Hamiltonian terms is the actual
  local term count in the non-TS model.
  - TS Trotter with resum(H) turns a compact 3-term TS Hamiltonian into 3N
  translated gates.
  - The optimized single-Pauli gate kernel is cheap, but applying O(N) gates per
  step plus repeated truncation dominates.
  - resum(H) materialization is minor; the expensive part is using the
  translated single-gate product formula as the computational primitive.

  So the issue is not merely expansion as storage. It is expansion as
  algorithmic work granularity.

  Implication for design

  The lazy resum change is insufficient and maybe not worth keeping as-is,
  because it removes allocation but not the algorithmic bottleneck.

  A useful generic next direction should target one of:

  1. Adaptive/coarser truncate_every, probably defaulting to something like 4–16
  or auto-tuned based on growth.
  2. A TS-native Trotter primitive that does not execute N single-site 
  translated gates as separate high-level operations.
  3. An anchored/internal representation, if we want to still represent ordered
  translated gates but reduce canonicalization/sorting overhead.
  4. Orbit-level product formula, but that starts becoming a different
  integrator unless done carefully.

  The immediate concrete low-risk performance lever is truncate_every; it gives
  a clear speedup in profiling without changing the mathematical gate order,
  only when truncation is applied during the step. But the deeper reason Trotter
  loses to RK4 is the O(N) translated-gate sweep, not resum(H) materialization.

✻ Sautéed for 18m 0s

❯ the fix suggested by you seems to be a bandaid fix rather than genuinely      
  exploiting the TS structure of the problem. Presumably, a TS native Trotter   
  primitive is a better avenue, which can lead to genuine gains. However, we    
  want to be sure that your hypothesis is correct. Formulate a definitive check 
  that verifies that this truely is the cause of the slow down. For now,        
  revert the changes made so we can pinpoint the exact problem. However, store  
  it somewhere so that we can pull this in when we are worried about finer      
  optimizations.                                                                

  Ran 2 shell commands (ctrl+o to expand)

● Reverted the working tree and stored the lazy prototype in a stash:

  stash@{0}: On 78: lazy-ts-trotter-prototype-before-ts-native

  Working tree is clean again. We can recover it later with:

  git stash apply stash@{0}

  or inspect it with:

  git stash show -p stash@{0}

  Definitive check for the real slowdown hypothesis

  The hypothesis we now need to test is:

  ▎ TS Trotter is slow because it uses a non-TS product formula over O(N) 
  ▎ translated single-Pauli gates. The issue is not resum(H) allocation, but 
  ▎ that the computational primitive has the wrong granularity: 6N-1 gate 
  ▎ updates per second-order step instead of a small number of TS-native 
  ▎ orbit/layer updates.

  A definitive check should isolate gate-count / translated-sweep granularity
  from all other effects.

  Check 1: scaling with system size at fixed TS complexity

  Use a Hamiltonian with fixed TS representative length, e.g. MFIM:

  length(Hts) = 3

  and choose truncation so the TS state size stays approximately comparable
  across N.

  Then measure one-step or short-run cost of:

  1. RK4 TS step:

  rk4(Hts, Ots, dt; truncation=...)
  2. current TS Trotter internal step:

  Hr = resum(Hts)
  gates = trotterize(Hr, dt; order=2)
  trotter_step!(representative(Ots), gates; truncation=...)

  Record:

  N
  length(Hts)
  length(resum(Hts))
  length(gates)
  length(Ots before)
  time RK4
  time Trotter
  time per Trotter gate
  trim calls

  Expected if hypothesis is true:

  length(Hts)         constant: 3
  length(gates)       grows ~ 6N
  RK4 cost            grows mildly/with TS state size
  Trotter cost        grows approximately linearly with gate count, all else
  equal
  Trotter/RK4 ratio   worsens with N

  This would show the bottleneck is not materialization but the O(N) gate sweep.

  Check 2: synthetic gate-list scaling at fixed state

  This is cleaner and more definitive.

  Take a fixed ordinary representative state Or of controlled size. Then apply
  only prefixes of the full Trotter gate list:

  gates = trotterize(resum(Hts), dt; order=2)

  for K in [16, 32, 64, 128, 256, length(gates)]
      Otest = copy(Or)
      trotter_step!(Otest, gates[1:K]; truncation=truncation)
      measure time
  end

  Also record:

  K
  time
  trim calls
  max intermediate length

  Expected if hypothesis is true:

  time ≈ linear/superlinear in K

  especially once branching and trimming are included.

  This directly tests whether the number of translated gates is the driver.

  Check 3: compare against an artificial “orbit primitive” cost floor

  We can mock the ideal TS-native primitive without implementing its physics.

  For each TS representative Hamiltonian term, do one TS commutator-like pass
  over relative shifts, similar to the existing TS commutator kernel, just to
  estimate the lower-bound cost of an orbit-level primitive.

  For MFIM, a second-order TS-orbit Strang step would have approximately:

  2 * length(Hts) - 1 = 5

  orbit factors instead of:

  2 * length(resum(Hts)) - 1 = 6N - 1

  So the computational-granularity ratio is:

  (6N - 1) / 5

  For N=64:

  383 / 5 ≈ 76.6x

  A real TS-native primitive would not be a single cheap Pauli rotation, so it
  will not get a full 76x improvement. But if profiling shows current Trotter
  time is dominated by hundreds of repeated gate sweeps/trims, then even a
  moderately more expensive orbit primitive can win.

  The definitive check here is not numerical correctness, but cost model
  validation:

  current Trotter = hundreds of ordinary gate passes
  RK4 TS = few TS-kernel passes

  If the time correlates with pass count, a TS-native primitive is the right
  avenue.

  Check 4: disable materialization cost entirely

  This was already partially done.

  Lazy gate generation reduced gate construction:

  N=64 lazy_gates_time ≈ 0.006s
  baseline resum + gates ≈ 0.16s one-time

  but total loop time remained ~seconds. Therefore materialization is not the
  main bottleneck.

  A definitive version is:

  Hr = resum(H)
  gates = trotterize(Hr, dt)
  # warm everything
  @time many_steps_using_cached_gates

  If cached-gate Trotter is still slow, resum is not the issue. We observed
  exactly that.

  The decisive benchmark I would run next

  I’d write one benchmark script that prints this table:

  N | len(Hts) | len(Hr) | gates | len(Ots) | RK4 step | Trotter step cached |
  Trotter trim time | gate-update time

  for e.g.

  N = 16, 32, 64, 128

  using the same M or adjusted M to keep TS state size comparable.

  Then we should see whether:

  Trotter step time / RK4 step time ≈ tracks gates / TS terms

  If yes, the cause is definitively the expanded translated-gate primitive, and
  the correct design is TS-native Trotter/orbit evolution rather than lazy
  generation or fold-frequency tweaks.

✻ Cogitated for 56s

❯ seems good. Try it out                                                        

● Write(~/.claude/jobs/513a65b3/tmp/definitive_ts_trotter_check.jl)
  ⎿  Wrote 168 lines to
     ../../../.claude/jobs/513a65b3/tmp/definitive_ts_trotter_check.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N, h)
       5     H = Operator(N)
       6     H += -h, "X", 1
       7     H += -h/2, "Z", 2
       8     H += "Z", 1, "Z", 2
       9     return OperatorTS{(N,)}(H)
      10 end
     … +186 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/definitive_ts_trotter_check.jl)
  ⎿  Added 1 line, removed 1 line
      128              trotter_step!(Ow, prefix; truncation=trunc)
      129          end
      130          Ow = copy(Or)
      131 -        trotter_step!(Ow, prefix; truncation=trim)            
      131 +        trotter_step!(Ow, prefix; truncation=o -> trim(o, max_
          +strings))                                                     
      132          @printf("K=%d time=%.6g perGate=%.6g trimCalls=%d trim
           Time=%.6g maxIn=%d finalLen=%d\n",
      133                  K, t, t/K, stats.calls, stats.time, stats.max_
           in, length(Ow))
      134          flush(stdout)

● Update(~/.claude/jobs/513a65b3/tmp/definitive_ts_trotter_check.jl)
  ⎿  Added 1 line, removed 1 line
      155              trotter_step!(Ow, gates; truncation=trunc, truncat
           e_every=te)
      156          end
      157          Ow = copy(Or)
      158 -        trotter_step!(Ow, gates; truncation=trim, truncate_eve
          -ry=te)                                                        
      158 +        trotter_step!(Ow, gates; truncation=o -> trim(o, max_s
          +trings), truncate_every=te)                                   
      159          @printf("te=%d time=%.6g trimCalls=%d trimTime=%.6g ma
           xIn=%d finalLen=%d\n",
      160                  te, t, stats.calls, stats.time, stats.max_in, 
           length(Ow))
      161          flush(stdout)

  Ran 1 shell command (ctrl+o to expand)

● I ran the checks on the reverted baseline. The results support the corrected
  hypothesis:

  ▎ The slowdown is mainly the O(N) translated single-gate sweep, not resum(H) 
  ▎ materialization. Current TS Trotter is using the non-TS single-Pauli 
  ▎ rotation primitive at the wrong granularity.

  1. Size-scaling check at fixed TS complexity

  Setup:

  - MFIM TS Hamiltonian, so length(Hts) = 3 for every N.
  - Grew Ots to about 2048 TS representative strings for each N.
  - trim cap 4096.
  - Compared one RK4 step against one cached current-style Trotter step using
  trotterize(resum(Hts), dt).

  Results:

  ┌─────┬──────────┬───────────────┬──────┬──────────┬───────┬───────┬─────┐
  │     │ length(H │ length(resum( │ Stra │ length(O │ RK4   │ Trott │ Rat │
  │  N  │   ts)    │     Hts))     │ ng g │   ts)    │ step  │  er   │ io  │
  │     │          │               │ ates │          │       │ step  │     │
  ├─────┼──────────┼───────────────┼──────┼──────────┼───────┼───────┼─────┤
  │  16 │        3 │            48 │   95 │     2048 │ 0.010 │ 0.017 │ 1.6 │
  │     │          │               │      │          │    5s │    2s │  4x │
  ├─────┼──────────┼───────────────┼──────┼──────────┼───────┼───────┼─────┤
  │  24 │        3 │            72 │  143 │     2048 │ 0.012 │ 0.019 │ 1.5 │
  │     │          │               │      │          │    8s │    2s │  1x │
  ├─────┼──────────┼───────────────┼──────┼──────────┼───────┼───────┼─────┤
  │  32 │        3 │            96 │  191 │     2048 │ 0.014 │ 0.031 │ 2.1 │
  │     │          │               │      │          │    7s │    9s │  7x │
  ├─────┼──────────┼───────────────┼──────┼──────────┼───────┼───────┼─────┤
  │  48 │        3 │           144 │  287 │     2048 │ 0.019 │ 0.032 │ 1.6 │
  │     │          │               │      │          │    3s │    4s │  8x │
  ├─────┼──────────┼───────────────┼──────┼──────────┼───────┼───────┼─────┤
  │  64 │        3 │           192 │  383 │     2048 │ 0.022 │ 0.047 │ 2.0 │
  │     │          │               │      │          │    6s │    0s │  8x │
  └─────┴──────────┴───────────────┴──────┴──────────┴───────┴───────┴─────┘

  One-time costs were tiny:

  N=64 resum(H)      ≈ 0.00027s
  N=64 gate building ≈ 0.000026s
  N=64 Trotter step  ≈ 0.047s

  So resum(H) allocation is absolutely not the main cost here.

  The Trotter cost is tied to applying 6N-1 gates. The scaling is not perfectly
  monotonic because the string content/commutation pattern changes with N, but
  the cost remains clearly controlled by the expanded translated-gate count.

  2. Prefix gate-list scaling check

  This is the most definitive test.

  Setup:

  - Fixed N=64
  - Fixed TS state size: length(Ots)=4096
  - Built full ordinary Strang gate list from resum(Hts): 383 gates
  - Applied only prefixes of length K

  Results:

  ┌──────────┬──────────┬───────────┬──────────────┐
  │ Prefix K │   Time   │ Time/gate │ Final length │
  ├──────────┼──────────┼───────────┼──────────────┤
  │        8 │ 0.00103s │ 0.000129s │         4704 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │       16 │ 0.00189s │ 0.000118s │         5039 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │       32 │ 0.00367s │ 0.000115s │         5724 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │       64 │ 0.00887s │ 0.000139s │         7548 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │      128 │  0.0430s │ 0.000336s │         8192 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │      256 │  0.0766s │ 0.000299s │         8192 │
  ├──────────┼──────────┼───────────┼──────────────┤
  │      383 │   0.107s │ 0.000281s │         8192 │
  └──────────┴──────────┴───────────┴──────────────┘

  This directly confirms the gate-sweep hypothesis.

  For small K, runtime scales almost linearly with gate count. Around K≈128, the
  operator hits the truncation cap and sorting/trimming effects kick in, but
  the dominant axis is still how many translated gates are applied.

  This is the exact computational primitive problem:

  TS Hamiltonian terms: 3
  current Trotter gates at N=64: 383

  So current TS Trotter is doing hundreds of non-TS single-gate passes to
  represent what should conceptually be a few TS orbit-level operations.

  3. truncate_every check

  Setup:

  - N=64
  - length(Ots)=4096
  - full 383-gate step
  - cap 8192

  Results:

  ┌────────────────┬─────────┬──────────────┬─────────────────────────┐
  │ truncate_every │  Time   │ Final length │ Max intermediate length │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │              1 │  0.108s │         8192 │                   12169 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │              2 │  0.102s │         8192 │                   13135 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │              4 │ 0.0958s │         8192 │                   13125 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │              8 │ 0.0926s │         8790 │                   19901 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │             16 │ 0.0918s │         9636 │                   21798 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │             32 │  0.103s │         9954 │                   32107 │
  ├────────────────┼─────────┼──────────────┼─────────────────────────┤
  │             64 │  0.191s │        14497 │                  135942 │
  └────────────────┴─────────┴──────────────┴─────────────────────────┘

  This shows truncate_every is indeed a tuning knob, but only a modest one. Best
  case here is roughly:

  0.108s -> 0.092s

  That is a ~15% improvement, not a structural fix.

  And if truncation is too infrequent, intermediate growth becomes bad. So this
  is a band-aid, not the core solution.

  Conclusion

  The definitive check validates:

  current TS Trotter slowdown = expanded translated-gate work granularity

  not:

  one-time resum(H) materialization

  and not primarily:

  folding/final contraction cost

  The problematic code path is exactly:

  - src/evolution.jl: TS Trotter uses Hr = resum(H)
  - src/trotter.jl: trotterize(Hr) creates 6N-1 Strang gates
  - src/trotter.jl: trotter_step! applies every gate sequentially, looping over
  the whole operator each time

  For N=64, RK4 performs a few TS-kernel passes; Trotter performs 383 ordinary
  gate passes.

  So yes: a genuine TS-native Trotter primitive is the right avenue. The
  primitive needs to operate at the level of TS Hamiltonian
  representatives/orbits, not individual translated ordinary Pauli rotations.
  The goal is to replace:

  6N - 1 ordinary gate sweeps

  with something closer to:

  2 * length(Hts) - 1 TS orbit sweeps

  for a second-order method, while preserving correctness of the finite-step
  product formula or defining a justified TS-preserving product formula.

✻ Baked for 2m 49s

❯ very good. Explicitly write down your findings in a markdown file. include    
  what you changed, and what results you got, and how this verifies your        
  hypothesis. it should be clear, but not verbose. This is meant only to        
  document our findings informally, so we can come back and verify this.        

  Searched for 1 pattern (ctrl+o to expand)

● Write(benchmark/ts_trotter_slowdown_findings.md)
  ⎿  Wrote 185 lines to benchmark/ts_trotter_slowdown_findings.md
      1 # TS Trotter slowdown notes
      2 
      3 Date: 2026-06-04
      4 
      5 ## Problem
      7 The current translation-symmetric (TS) Trotter path does not beat
         RK4 on the TS mixed-field Ising example, unlike the non-TS case.
      8 
      9 Initial suspicion was that the slowdown came mainly from expandin
        g `H` with `resum(H)` and evolving outside the TS representation.
     … +199 lines (ctrl+o to expand)

● Wrote the findings to:

  benchmark/ts_trotter_slowdown_findings.md

  It includes:

  - current baseline code path,
  - what the lazy prototype changed,
  - where the prototype is stored,
  - before/after timings,
  - component profile results,
  - definitive scaling/prefix/truncation checks,
  - and the conclusion that the real issue is the O(N) translated single-gate
  primitive, not resum(H) materialization.

✻ Brewed for 53s

※ recap: We’re investigating why TS Trotter is slower than RK4 and have documented evidence that the bottleneck is the expanded translated-gate primitive. Next, design a TS-native Trotter/orbit primitive and validate it against the benchmark. (disable recaps in /config)
❯ hi                                                                            

● Hi! Want to continue from the TS-native Trotter design, or switch to something
  else?

✻ Cooked for 4s

※ recap: We’re investigating why TS Trotter is slower than RK4 and documented that the bottleneck is the expanded translated-gate primitive. Next, design a genuinely TS-native Trotter/orbit primitive instead of optimizing lazy resum.

❯ /plan find different avenues for a TS native trotter primitive                
  ⎿  Enabled plan mode

● Plan(Design TS Trotter avenues) Sonnet 4.6
  ⎿  Done (0 tool uses · 0 tokens · 0s)
  (ctrl+o to expand)

● Plan(Design TS Trotter avenues) Haiku 4.5
  ⎿  Backgrounded agent (↓ to manage · ctrl+o to expand)

● Agent "Design TS Trotter avenues" completed · 0s

● I’ll map the design space first and keep it generic-first, then propose which
  avenue is most worth implementing. I’ve also launched a separate read-only
  design pass to sanity-check the options before I finalize the plan.

✻ Sautéed for 51s

❯ find different avenues for a TS native trotter primitive                      

● Task Output a34b4a509b1e9c73f
  ⎿  Error: No task found with ID: a34b4a509b1e9c73f

● Yes — if we want genuine gains, the primitive should stop treating a TS
  Hamiltonian as O(N) ordinary gates. Below are the main avenues I see for a
  TS-native Trotter primitive, ordered roughly from most generic to more
  structure-exploiting.

  I’ll use:

  H = \sum_{a=1}^m H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_{r \in G} T_r P_a T_r^\dagger

  where each P_a is one representative TS Hamiltonian term.

  A TS-native second-order Trotter should ideally split over the orbit 
  Hamiltonians:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2} L_{H_1^{orb}}}
  e^{\frac{\Delta t}{2} L_{H_2^{orb}}}
  \cdots
  e^{\Delta t L_{H_m^{orb}}}
  \cdots
  e^{\frac{\Delta t}{2} L_{H_2^{orb}}}
  e^{\frac{\Delta t}{2} L_{H_1^{orb}}},

  where:

  L_H(O) = i[H, O] / \hbar.

  This has only 2m-1 TS-preserving factors, instead of 2 * length(resum(H)) - 1 
  = O(Nm) ordinary gates.

  For MFIM, this is the difference between:

  5 orbit factors

  and:

  383 ordinary gates at N=64

  ---
  Avenue 1: Matrix-free orbit Liouvillian exponential

  Idea

  For each orbit term H_a^{orb}, define a TS-native linear operator:

  L_a(O) = i[H_a^{orb}, O] / \hbar.

  Then apply:

  O \mapsto e^{\tau L_a} O

  directly in OperatorTS using a matrix-free exponential method: Krylov/Arnoldi,
  Lanczos-like, Chebyshev, etc.

  The outer Trotter is still second-order Strang over orbit Hamiltonians.

  Code fit

  Introduce something like:

  struct OrbitLiouvillian{...}
      term::PauliStringTS
      coeff
      hbar
  end

  apply(L::OrbitLiouvillian, O::OperatorTS)

  Internally this can initially call the existing TS commutator path:

  commutator(OperatorTS(term), O)

  Then later optimize this with support-aware relative-shift enumeration.

  Then:

  orbit_expv(L, O, τ; tol, krylov_dim)

  computes exp(τL)O.

  Cost

  If Krylov dimension is q, cost per Trotter step is roughly:

  (2m - 1) * q * cost(one TS orbit commutator)

  For MFIM, m = 3, so:

  5q TS commutator-like passes

  RK4 costs roughly:

  4 full-H TS commutator passes

  where each full-H pass loops over all m terms.

  So this only wins if:

  - q is modest,
  - orbit commutator application is optimized,
  - larger dt is possible because orbit exponentials are handled more accurately
  than RK4 stages.

  Correctness

  This is clean mathematically:

  - Strang over H_a^{orb} gives second-order Trotter error.
  - Inner expv tolerance controls the numerical error of each orbit exponential.
  - Every substep preserves TS exactly, up to truncation.

  Risks

  - If q is large, this becomes RK-like or slower.
  - Needs a good matrix-free exponential implementation.
  - Requires clear truncation policy inside Krylov/exponential action.

  Validation

  - Small N: compare against dense exact evolution.
  - Check second-order scaling: global error should scale like dt^2.
  - Check first derivative:

  (e^{\epsilon L_a}O - O)/\epsilon \approx L_a O.
  - Benchmark cost versus current Trotter and RK4 for increasing N.

  This is probably the most generic TS-native avenue.

  ---
  Avenue 2: Exact commuting-orbit primitive

  Idea

  Some orbit Hamiltonians have mutually commuting translated copies:

  [P_a(r), P_a(r')] = 0
  \quad \forall r,r'.

  Then:

  e^{i \tau H_a^{orb}}
  =
  \prod_r e^{i \tau c_a P_a(r)}.

  Current code applies this as O(N) ordinary gates. But if the orbit is
  commuting, we can apply the whole orbit action directly to a TS operator.

  For each TS string Q, only translated copies of P_a that anticommute with Q
  matter. For local strings, the number of such relative shifts is usually
  small, independent of total N.

  So instead of looping over all translations:

  for r in 1:N
      apply single gate P_a(r)
  end

  we do:

  find active relative shifts that anticommute with Q
  apply combined commuting rotation
  fold result into TS representatives

  Code fit

  Add generic detection:

  is_commuting_orbit(p::PauliStringTS)

  Then:

  orbit_rotation_commuting(p_ts, coeff, O_ts, τ)

  This would be an exact primitive for that orbit.

  Cost

  Potentially:

  O(length(Ots) * active_relative_shifts)

  instead of:

  O(N * length(Oordinary))

  For local Hamiltonian terms and truncated strings, active_relative_shifts
  should be closer to local support size than system size.

  Correctness

  Exact for each commuting orbit.

  The outer Strang formula over orbit Hamiltonians remains second order.

  Risks

  - Only applies when translated copies commute.
  - Applying all active commuting rotations can branch as 2^k, where k is number
  of anticommuting relative shifts. For local/truncated strings this may be
  small; for large support strings it may explode.
  - Needs careful coefficient/phase handling.

  Relevance

  MFIM terms are favorable:

  - onsite X orbit: commuting,
  - onsite Z orbit: commuting,
  - nearest-neighbor ZZ orbit: commuting.

  So this would likely give the benchmark a real speedup.

  But it is still generic in the sense that it detects the property; it is not
  hard-coded to MFIM.

  ---
  Avenue 3: Generic graph-colored orbit primitive

  Idea

  For a general orbit term, translated copies may not all commute. Build the
  conflict graph:

  vertices = translations r
  edge r--s if [P_a(r), P_a(s)] != 0

  Then color the graph. Each color class contains mutually commuting translated
  copies.

  Approximate:

  e^{\tau L_{H_a^{orb}}}

  by a Strang product over color classes:

  e^{\frac{\tau}{2} L_{C_1}}
  e^{\frac{\tau}{2} L_{C_2}}
  \cdots
  e^{\tau L_{C_k}}
  \cdots
  e^{\frac{\tau}{2} L_{C_2}}
  e^{\frac{\tau}{2} L_{C_1}}.

  Each C_j is a commuting suborbit/layer and can use the commuting-orbit
  primitive above.

  Code fit

  orbit_coloring(p_ts) -> Vector{ColorClass}

  Then:

  orbit_exp_colored(p_ts, coeff, O_ts, τ)

  Internally applies commuting layers.

  Cost

  If the number of colors is bounded by local support, cost becomes roughly:

  O(num_colors * cost(commuting orbit layer))

  rather than O(N) gates.

  For local 1D nearest-neighbor terms this recovers even/odd-like behavior, but
  generated generically.

  Correctness

  This adds an inner Trotter approximation for each noncommuting orbit
  exponential.

  The full method becomes nested Strang:

  1. outer Strang over Hamiltonian representative orbits,
  2. inner Strang over color classes within each orbit.

  It can still be second order if composed symmetrically.

  Risks

  - Individual color classes may not preserve full one-site translation
  symmetry. They may preserve only a subgroup, e.g. even-bond layers.
  - Plain OperatorTS may be insufficient for intermediate states.
  - May require either:
    - sublattice-TS representations,
    - an anchored TS workspace,
    - or temporary ordinary representation.

  So this is promising but more complex than the commuting-orbit case.

  ---
  Avenue 4: Anchored TS workspace for ordered products

  Idea

  Current OperatorTS quotients out absolute translation origin. That is why it
  cannot represent intermediate states after one translated gate: a single
  translated gate breaks translation symmetry.

  Introduce an internal representation:

  (canonical TS representative, anchor/origin/frame) -> coefficient

  This keeps enough origin information to apply ordered translated gates without
  fully expanding into ordinary strings.

  Then a Trotter step becomes:

  OperatorTS
    -> AnchoredTSOperator
    -> apply ordered translated gates in quotient form
    -> fold anchors away
    -> OperatorTS

  Code fit

  New internal type only; public API still returns OperatorTS.

  struct AnchoredTSOperator
      strings::Vector{PauliStringTS}
      anchors::Vector{Shift}
      coeffs::Vector
  end

  Cost

  This probably reduces duplication and canonicalization overhead, but it does
  not remove the O(N) translated-gate sweep by itself.

  So it is not the main structural fix, but it may be necessary for
  graph-colored or symmetrized methods.

  Correctness

  Can preserve the exact current product formula if implemented carefully.

  Risks

  - More complicated representation and hashing.
  - May only provide constant-factor improvements.
  - Does not by itself solve the 6N-1 gate problem.

  I would treat this as infrastructure, not the primary primitive.

  ---
  Avenue 5: Cached TS superoperator for a fixed truncated basis

  Idea

  For a fixed active TS basis B, build the matrix of a whole TS-native step:

  S_{TS}: B \to B.

  Then each time step is just:

  coeffs_next = S_TS * coeffs

  This could be built for either:

  - the orbit-Liouvillian exponential,
  - the full Strang step,
  - or an approximate projected step.

  Code fit

  This is more like a separate engine:

  build_ts_step_operator(Hts, basis, dt; method)
  apply_step(S, coeffs)

  Cost

  High precompute cost, low per-step cost.

  Good when:

  - the basis is fixed,
  - many time steps are needed,
  - many observables/runs reuse the same Hamiltonian.

  Bad when:

  - trim changes the active basis adaptively,
  - basis size is large,
  - one-off evolution.

  Correctness

  Depends on the underlying step used to build the superoperator.

  Projection/truncation must be explicit.

  Risks

  - Memory can blow up.
  - Adaptive truncation conflicts with a fixed matrix.
  - Less idiomatic for the current dynamic Pauli-string representation.

  Probably not the first implementation.

  ---
  Avenue 6: Group-averaged Trotter step

  Idea

  Given an ordinary ordered product formula S, define a TS-preserving map:

  \Phi(O)
  =
  \frac{1}{|G|}
  \sum_{g \in G}
  T_g^\dagger S(T_g O T_g^\dagger) T_g.

  For TS input O, this is equivalent to averaging all translated choices of
  product-formula origin/order.

  This produces an explicitly TS-preserving finite-step map.

  Code fit

  Could be represented as:

  symmetrized_trotter_step(Hts, Ots, dt)

  Cost

  Naively terrible:

  O(N) full Trotter sweeps

  Sampling could reduce cost, but then the method becomes stochastic.

  Correctness

  The averaged map is TS-preserving and inherits the product formula order if
  the base formula has that order.

  Risks

  - Too expensive exactly.
  - Stochastic version changes semantics and introduces variance.
  - Not likely to be the performance solution.

  This is useful conceptually, but probably not practical.

  ---
  Avenue 7: Support-aware TS commutator/orbit kernel

  This is not a full primitive by itself, but it is probably essential for any
  TS-native primitive.

  Idea

  The current TS binary kernel loops over all translations:

  for s in all_shifts(...)
      p, k = op(rep1, shift(rep2, ..., s))
  end

  But for local Pauli strings, most shifts commute because supports do not
  overlap.

  We can enumerate only relative shifts where supports overlap or can
  anticommute.

  Code fit

  Add helper:

  active_relative_shifts(p_rep, q_rep, Ls, Ps)

  Then use it in:

  - TS commutator,
  - orbit Liouvillian apply,
  - commuting-orbit primitive,
  - graph construction.

  Cost

  Could reduce TS kernel cost from:

  O(N)

  relative shifts per pair to roughly:

  O(local support overlap candidates)

  for local/truncated operators.

  Correctness

  Exact if active shifts include every possible nonzero commutator/product
  contribution.

  Risks

  - More involved for multidimensional lattices and periodic/open flags.
  - For dense/high-weight strings, active shifts may still be large.
  - Need careful treatment of stabilizers/periodic representatives.

  This is likely a prerequisite optimization for the best TS-native approaches.

  ---
  My ranking

  Most promising generic path

  1. Orbit-Liouvillian exponential + support-aware TS kernel

  This is the cleanest genuinely TS-native primitive:

  Strang over TS orbit Hamiltonians
  each orbit exponential applied in OperatorTS

  It is generic and mathematically clean.

  Main question: can orbit exponential action be made cheaper than RK4? That
  depends on Krylov dimension and kernel speed.

  Most promising benchmark-winning path

  2. Exact commuting-orbit primitive with generic detection

  This should directly attack MFIM without hard-coding MFIM.

  It replaces hundreds of translated gates by a few exact commuting orbit
  rotations.

  This is less universal than Avenue 1, but still generic because it applies
  whenever the orbit has commuting translated copies.

  Most useful fallback/infrastructure

  3. Graph-colored orbit primitive + anchored/subgroup TS workspace

  This handles noncommuting translated copies more Trotter-like, but
  implementation complexity is higher.

  ---
  Recommended exploration order

  1. Prototype support-aware orbit commutator application.
  Verify it matches current commutator(OperatorTS(single_term), Ots) exactly but
  avoids scanning all shifts where possible.
  2. Prototype exact commuting-orbit primitive.
  Use generic commutation detection. Validate on MFIM terms. This is likely the
  quickest route to a real TS Trotter speedup.
  3. Prototype orbit-Liouvillian expv.
  Use the orbit commutator as a matrix-free operator. Compare Krylov dimensions
  and errors against dense exact for small N.
  4. Only then consider graph coloring / anchored TS.
  Use this if important Hamiltonian terms have noncommuting orbits and Krylov is
  too slow.

  The key design shift is:

  do not Trotterize over translated ordinary Pauli strings

  Instead Trotterize over:

  translation orbits of Pauli strings

  and make the primitive:

  exp(τ L_{orbit}) acting on OperatorTS

  rather than:

  for every translation: apply one ordinary Pauli rotation

✻ Sautéed for 2m 52s

❯ while I understand that MFIM may benefit from commuting orbit primitives, we  
  will not focus on exploiting the structure of MFIM. we will only exploit the  
  TS structure. This is the goal of this exercise. "Design a second order       
  Trotter method that exploits translation invariance and outperform RK4 in the 
  TS case." this is our prompt                                                  

● Understood. We should not target MFIM’s commuting structure. The design must
  exploit only translation invariance, i.e. the fact that the Hamiltonian is
  represented as a small number of translation orbits.

  So the core primitive should be:

  H = \sum_{a=1}^m H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_r T_r P_a T_r^\dagger

  and the second-order Trotter method should split over the TS orbit 
  Hamiltonians:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2}L_{H_1^{orb}}}
  e^{\frac{\Delta t}{2}L_{H_2^{orb}}}
  \cdots
  e^{\Delta t L_{H_m^{orb}}}
  \cdots
  e^{\frac{\Delta t}{2}L_{H_2^{orb}}}
  e^{\frac{\Delta t}{2}L_{H_1^{orb}}},

  where:

  L_H(O) = i[H, O] / \hbar.

  This uses only TS structure: the Hamiltonian is split into translation orbits,
  not into individual translated Pauli strings.

  The hard part is applying one orbit flow:

  O \mapsto e^{\tau L_{H_a^{orb}}} O

  without expanding into all translated single-Pauli rotations.

  Below are the TS-only avenues worth considering.

  ---
  Avenue 1: Matrix-free orbit exponential in OperatorTS

  This is the cleanest TS-native design.

  Primitive

  Define a matrix-free orbit Liouvillian:

  L_a(Ots) = im / hbar * commutator(Ha_orbit, Ots)

  where:

  Ha_orbit = OperatorTS(P_a) * c_a

  Then apply:

  exp(τ * L_a) * Ots

  directly in OperatorTS, using Krylov/Arnoldi or another exponential-action
  method.

  Why this is TS-native

  Every operation is in the TS algebra:

  OperatorTS -> OperatorTS

  No translated ordinary gates. No resum(H). No dependence on whether translated
  copies commute.

  Expected cost

  Outer Strang uses:

  2m - 1 orbit exponentials

  For each orbit exponential, Krylov dimension q gives about:

  q orbit-commutator applications

  So cost per step is roughly:

  (2m - 1) * q orbit-commutators

  RK4 costs:

  4 full-H commutators = 4m orbit-term passes

  For MFIM, m = 3:

  RK4:          ~12 orbit-term passes
  Orbit Strang: 5q orbit-term passes

  So this wins only if q is small, e.g. q = 1–2, or if the orbit commutator
  kernel is faster/sparser than the full-H RK4 intermediates.

  Correctness

  Clean:

  - Strang over orbit Hamiltonians is second order.
  - Krylov/exponential tolerance controls inner error.
  - Translation symmetry is preserved at every substep.

  Main risk

  If Krylov dimension is more like q=4–8, this will not beat RK4.

  Definitive check before full implementation

  Implement only a prototype orbit expv wrapper around the existing TS
  commutator and measure:

  q needed for small-N dense accuracy
  cost per orbit expv
  cost of one full orbit-Strang step

  If q <= 2 is enough at useful dt, this avenue is promising. If not, probably
  not.

  ---
  Avenue 2: Second-order polynomial orbit exponential

  This is a cheaper, more approximate TS-native primitive.

  Primitive

  Approximate each orbit exponential by the second-order Taylor action:

  e^{\tau L_a} O
  \approx
  O + \tau L_a O + \frac{\tau^2}{2} L_a^2 O.

  Then use the same outer Strang ordering over orbit Hamiltonians.

  Why this is TS-native

  Again:

  OperatorTS -> OperatorTS

  Only uses TS commutators.

  No translated gates, no MFIM structure.

  Expected cost

  Each orbit factor costs two orbit-commutator applications:

  L_a O
  L_a(L_a O)

  So one second-order orbit-Strang step costs:

  2 * (2m - 1) orbit-commutators

  For MFIM:

  2 * 5 = 10 orbit-commutator applications

  RK4 costs about:

  4 * 3 = 12 orbit-term commutator applications

  So this is in the right cost range to plausibly beat RK4, especially if
  intermediates remain smaller.

  Correctness

  This should be globally second order if the orbit exponential approximation is
  second-order accurate and used consistently inside the symmetric split.

  But it is no longer an exact product of exponentials. It is more like:

  Strang splitting over orbit flows, with second-order TS-native subflow
  approximations

  So we should be careful about calling it “Trotter” without qualification.

  Pros

  - Simple.
  - Fully TS-native.
  - No Krylov implementation initially.
  - Directly tests whether orbit-level splitting helps.

  Cons

  - Does not preserve exact unitary conjugation per orbit.
  - May have worse stability than exact orbit exponential.
  - Accuracy could be closer to RK2 than Trotter if dt is too large.

  Validation

  - Dense small-N comparison.
  - Verify global error scales as O(dt^2).
  - Compare against RK4 at matched accuracy, not only matched dt.

  This is likely the fastest prototype to test the TS-orbit splitting idea.

  ---
  Avenue 3: Cayley / implicit midpoint orbit flow

  This is a more stable second-order TS-native orbit primitive.

  Primitive

  For each orbit Liouvillian L_a, use the Cayley transform:

  e^{\tau L_a}
  \approx
  \left(I - \frac{\tau}{2}L_a\right)^{-1}
  \left(I + \frac{\tau}{2}L_a\right).

  This is second-order and time-reversible.

  Why this is TS-native

  The only required operation is applying L_a in OperatorTS. The solve is
  performed in the TS coefficient space.

  Expected cost

  Depends on iterative solve iterations s.

  Cost per orbit factor:

  s orbit-commutator applications

  Total:

  (2m - 1) * s

  Pros

  - Better stability than Taylor.
  - Symmetric/time-reversible.
  - Natural for skew-Hermitian Liouvillians.

  Cons

  - Requires an iterative linear solver over adaptive sparse Pauli-string
  operators.
  - Truncation inside the solve is delicate.
  - Might be slower than RK4 unless iteration count is tiny.

  Validation

  Same as Krylov, plus check solver residual and norm behavior.

  This is probably not the first thing to implement, but it is a principled
  TS-native route.

  ---
  Avenue 4: TS projected exponential basis / adaptive local TS basis

  This is a more ambitious orbit-native approach.

  Primitive

  For the current active TS operator basis B, build a local matrix
  representation of L_a restricted to the truncated basis plus generated
  neighbors.

  Then apply:

  e^{\tau L_a}

  as a small sparse/dense exponential on that adaptive TS basis.

  Why this is TS-native

  The basis consists of PauliStringTS representatives only.

  Expected cost

  High setup cost, potentially low repeated step cost if the basis stabilizes.

  Pros

  - Direct orbit exponential.
  - Could be very fast for repeated many-step runs.
  - Natural place to reuse basis information between steps.

  Cons

  - Basis changes under trimming.
  - Building enough closure may be expensive.
  - Risk of dense matrices or large fill-in.

  This is probably useful only if Krylov/Taylor shows the orbit split is
  promising.

  ---
  Avenue 5: Symmetrized group-averaged Trotter map

  This is TS-only but likely expensive.

  Primitive

  Take an ordinary translated-gate product formula S, then average over
  translated origins:

  \Phi(O)
  =
  \frac{1}{|G|}
  \sum_{g \in G}
  T_g^\dagger S(T_g O T_g^\dagger) T_g.

  For TS input, this defines a TS-preserving map.

  Why this is TS-native conceptually

  The resulting map acts on the TS subspace and does not privilege one origin.

  Problem

  Naively this costs O(N) full Trotter sweeps, worse than the current method.

  It is useful as a mathematical reference, not a performance path.

  ---
  Important limitation

  Translation symmetry alone gives us a compact representation and efficient
  commutator algebra.

  It does not automatically give a closed-form two-branch Pauli rotation for:

  e^{\tau L_{H_a^{orb}}}

  because H_a^{orb} is generally a sum of noncommuting translated Pauli strings.

  The non-TS Trotter speedup comes from this identity for a single Pauli
  generator:

  e^{i\theta P/2} Q e^{-i\theta P/2}

  which maps one Pauli string to at most two Pauli strings.

  There is no equivalent generic two-branch rule for an entire translation orbit
  unless we exploit extra algebraic structure. Since we are not allowed to rely
  on such structure, the TS-native primitive must be an orbit Liouvillian 
  exponential/action method, not a direct single-Pauli rotation analogue.

  That is the core design constraint.

  ---
  Recommended plan

  Step 1: Define the target method

  Introduce conceptually:

  TrotterTSOrbit(; order=2, substep=:taylor2 or :krylov)

  with outer Strang over TS Hamiltonian representatives.

  For H::OperatorTS with m terms, the step is:

  half orbit 1
  half orbit 2
  ...
  full middle orbit
  ...
  half orbit 2
  half orbit 1

  Each orbit factor calls:

  orbit_flow(term, coeff, Ots, τ)

  and returns OperatorTS.

  Step 2: Prototype orbit_flow with Taylor-2

  This is the cheapest proof of concept.

  orbit_flow_a(O, τ)
  =
  O + τL_aO + τ^2/2 L_a^2O

  Then benchmark against RK4 and current Trotter.

  This tells us whether orbit-level TS splitting is actually promising before
  investing in Krylov or Cayley.

  Step 3: If Taylor-2 has the right cost profile, add Krylov expv

  Use the same L_a application but improve the subflow accuracy/stability.

  Compare:

  Taylor-2 orbit Strang
  Krylov orbit Strang
  RK4
  current expanded Trotter

  at matched error.

  Step 4: Optimize the orbit commutator kernel

  The existing TS commutator kernel loops over all translations. That still uses
  TS structure, but may be more expensive than necessary.

  A support-aware relative-shift kernel would help all TS-native orbit
  approaches.

  ---
  What I would not do as the main path

  - Do not focus on commuting-orbit exact rotations.
  - Do not focus on even/odd or graph coloring.
  - Do not optimize lazy resum(H) first.
  - Do not treat truncate_every tuning as the solution.

  Those may be later optimizations, but they do not address the prompt’s real
  target: a second-order method that exploits translation invariance
  generically.

  ---
  Bottom line

  The most aligned answer to the prompt is:

  ▎ A second-order Strang splitting over translation-orbit Hamiltonians, where 
  ▎ each orbit exponential is applied directly in OperatorTS using a TS-native 
  ▎ matrix-free orbit Liouvillian action.

  The first practical implementation should use the simplest TS-native orbit
  flow approximation, probably Taylor-2, to test whether orbit-level splitting
  beats RK4. If it does, upgrade the subflow to Krylov/Cayley for accuracy and
  stability.

✻ Baked for 1m 55s

❯ go ahead                                                                      

● I’ll implement the generic TS-orbit Strang prototype with a Taylor-2 orbit
  subflow first, then run correctness/error-scaling and timing checks.

  Read 3 files (ctrl+o to expand)

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line, removed 1 line
      24  export Circuits
      25  export pauli_weight, support, pauli_rotation
      26  export TrotterGate, trotterize, trotter_step!
      27 -export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
         -r, RK4, DOPRI5, Exact                                          
      27 +export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
         +r, TrotterTS, RK4, DOPRI5, Exact                               
      28  
      29  using Random
      30  using LinearAlgebra

● Update(src/evolution.jl)  ⎿  Added 15 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 +    TrotterTS(; order=2, substep=:taylor2)                     
      30 +                                                               
      31 +Translation-symmetric orbit-level product formula. For `H::Oper
         +atorTS`, splits the                                            
      32 +Hamiltonian into its representative translation orbits and appl
         +ies a first-order                                              
      33 +(`order=1`) or second-order (`order=2`, Strang) product over th
         +ose TS orbit                                                   
      34 +Liouvillians. The current `substep=:taylor2` applies each orbit
         + flow with a                                                   
      35 +second-order Taylor approximation in the TS operator algebra.  
      36 +"""                                                            
      37 +struct TrotterTS <: AbstractEvolutionMethod                    
      38 +    order::Int                                                 
      39 +    substep::Symbol                                            
      40 +end                                                            
      41 +TrotterTS(; order::Integer=2, substep::Symbol=:taylor2) = Trott
         +erTS(Int(order), substep)                                      
      42 +                                                               
      43 +"""                                                            
      44      RK4()
      45  
      46  Classical fixed-step 4th-order Runge–Kutta. Takes one internal 
          step per

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      126  
      127  # Keyword arguments
      128  - `method::AbstractEvolutionMethod = RK4()`. One of [`Trotter`
           ](@ref),
      129 -  [`RK4`](@ref), [`DOPRI5`](@ref), [`Exact`](@ref).           
      129 +  [`TrotterTS`](@ref), [`RK4`](@ref), [`DOPRI5`](@ref), [`Exac
          +t`](@ref).                                                    
      130  - `truncation`: function `O -> O` applied after every internal
            step. Default
      131    `identity`.
      132  - `dissipation`: function `(O, dt) -> O` applied after every i
           nternal step. The

● Update(src/evolution.jl)
  ⎿  Added 63 lines
      294      return EvolutionResult(collect(tspan), history, OperatorTS
           {Ls,Ps}(Or))
      295  end
      296  
      297 +function _orbit_liouvillian(Ha::Operator{<:PauliStringTS}, O::
          +Operator{<:PauliStringTS}, hbar::Real)                        
      298 +    return 1im * commutator(Ha, O) / hbar                     
      299 +end                                                           
      300 +                                                              
      301 +function _orbit_flow_taylor2(Ha::Operator{<:PauliStringTS}, O:
          +:Operator{<:PauliStringTS}, dt::Real, hbar::Real, truncation) 
      302 +    k1 = truncation(_orbit_liouvillian(Ha, O, hbar))          
      303 +    k2 = truncation(_orbit_liouvillian(Ha, k1, hbar))         
      304 +    return truncation(O + dt * k1 + (dt^2 / 2) * k2)          
      305 +end                                                           
      306 +                                                              
      307 +function _orbit_terms(H::Operator{<:PauliStringTS})           
      308 +    terms = typeof(H)[]                                       
      309 +    for j in 1:length(H)                                      
      310 +        c, p = H[j]                                           
      311 +        push!(terms, typeof(H)([p], [H.coeffs[j]]))           
      312 +    end                                                       
      313 +    return terms                                              
      314 +end                                                           
      315 +                                                              
      316 +function _trotterts_step(Hterms, O::Operator{<:PauliStringTS},
          + dt::Real, method::TrotterTS, hbar::Real, truncation)         
      317 +    method.order ∈ (1, 2) || throw(ArgumentError("order must b
          +e 1 or 2, got $(method.order)"))                              
      318 +    method.substep === :taylor2 || throw(ArgumentError("only s
          +ubstep=:taylor2 is currently implemented"))                   
      319 +    isempty(Hterms) && return O                               
      320 +    if method.order == 1 || length(Hterms) == 1               
      321 +        for Ha in Hterms                                      
      322 +            O = _orbit_flow_taylor2(Ha, O, dt, hbar, truncatio
          +n)                                                            
      323 +        end                                                   
      324 +    else                                                      
      325 +        L = length(Hterms)                                    
      326 +        for j in 1:(L - 1)                                    
      327 +            O = _orbit_flow_taylor2(Hterms[j], O, dt / 2, hbar
          +, truncation)                                                 
      328 +        end                                                   
      329 +        O = _orbit_flow_taylor2(Hterms[L], O, dt, hbar, trunca
          +tion)                                                         
      330 +        for j in (L - 1):-1:1                                 
      331 +            O = _orbit_flow_taylor2(Hterms[j], O, dt / 2, hbar
          +, truncation)                                                 
      332 +        end                                                   
      333 +    end                                                       
      334 +    return O                                                  
      335 +end                                                           
      336 +                                                              
      337 +function _evolve(method::TrotterTS, H::Operator{<:PauliStringT
          +S}, O::Operator{<:PauliStringTS}, tspan;                      
      338 +                 truncation, dissipation, fout, hbar)         
      339 +    qubitsize(H) == qubitsize(O) && periodicflags(H) == period
          +icflags(O) ||                                                 
      340 +        throw(DimensionMismatch("H and O must share the same t
          +ranslation-symmetry lattice"))                                
      341 +    n = length(tspan)                                         
      342 +    history = _alloc_history(fout, O, n)                      
      343 +    Hterms = _orbit_terms(H)                                  
      344 +    O = copy(O)                                               
      345 +    for i in ProgressBar(1:(n - 1))                           
      346 +        dt = tspan[i + 1] - tspan[i]                          
      347 +        O = _trotterts_step(Hterms, O, dt, method, hbar, trunc
          +ation)                                                        
      348 +        O = dissipation(O, dt)                                
      349 +        O = truncation(O)                                     
      350 +        _save!(history, fout, O, i + 1)                       
      351 +    end                                                       
      352 +    return EvolutionResult(collect(tspan), history, O)        
      353 +end                                                           
      354 +                                                              
      355 +function _evolve(::TrotterTS, H::AbstractOperator, O::Abstract
          +Operator, tspan;                                              
      356 +                 truncation, dissipation, fout, hbar)         
      357 +    throw(ArgumentError("TrotterTS evolution via `evolve` is i
          +mplemented for `OperatorTS` only, not for `$(typeof(H))`."))  
      358 +end                                                           
      359 +                                                              
      360  function _evolve(::Trotter, H::AbstractOperator, O::AbstractOp
           erator, tspan;
      361                   truncation, dissipation, fout, hbar)
      362      throw(ArgumentError("Trotter evolution via `evolve` is imp
           lemented for `Operator` and " *

  Ran 2 shell commands (ctrl+o to expand)

● Update(test/evolution.jl)
  ⎿  Added 3 lines
      32      res3 = evolve(H, O0, times; method = DOPRI5(), fout=fout, t
          runcation=truncation).history
      33      res4 = evolve(resum(H), resum(O0), times; method = Trotter(
          ), fout=fout_trotter, truncation=truncation).history
      34      res5 = evolve(H, O0, times; method = Trotter(), fout=fout, 
          truncation=truncation).history
      35 +    res6 = evolve(H, O0, times; method = TrotterTS(), fout=fout
         +, truncation=truncation).history                               
      36  
      37      @test norm(res2 .- res1)/norm(res1) < 1e-7
      38      @test norm(res3 .- res1)/norm(res1) < 1e-7
      39      @test norm(res4 .- res1)/norm(res1) < 1e-5
      40      @test norm(res5 .- res1)/norm(res1) < 1e-4
      41 +    @test norm(res6 .- res1)/norm(res1) < 1e-3                 
      42 +    @test evolve(H, O0, [0, dt]; method = TrotterTS(), fout=not
         +hing, truncation=truncation).final isa typeof(O0)              
      43  
      44      dissipation(O, dt) = add_noise(O, 0.1*dt)
      45      res1 = evolve(H, O0, times; method = RK4(), fout=fout, diss
          ipation=dissipation).history

  Ran 2 shell commands (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/compare_trotterts_trunc.jl)
  ⎿  Wrote 35 lines to
     ../../../.claude/jobs/513a65b3/tmp/compare_trotterts_trunc.jl
       1 using PauliStrings
       2 
       3 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       5 end
       6 function Xtot(N)
       7     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       8 end
      10     [typeof(H)([H.strings[j]],[H.coeffs[j]]) for j in 1:length(H
         )]
     … +33 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/compare_trotterts_trunc.jl)
  ⎿  Added 2 lines, removed 2 lines
      27      H=MFIM(N,0.5); O0=Xtot(N); terms=orbit_terms(H); trunc(o)=t
          rim(o,2^M); dt=0.1
      28      # warm
      29      step(terms,O0,dt,trunc,flow_stage); step(terms,O0,dt,trunc,
          flow_end)
      30 -    GC.gc(); t1=@elapsed begin O=O0; for _ in 1:steps; O=step(t
         -erms,O,dt,trunc,flow_stage); O=trunc(O); end; println((stage_le
         -n=length(O))); end                                             
      31 -    GC.gc(); t2=@elapsed begin O=O0; for _ in 1:steps; O=step(t
         -erms,O,dt,trunc,flow_end); O=trunc(O); end; println((end_len=le
         -ngth(O))); end                                                 
      30 +    GC.gc(); t1=@elapsed begin O=O0; for _ in 1:steps; O=step(t
         +erms,O,dt,trunc,flow_stage); O=trunc(O); end; println((stage_le
         +n=length(O),)); end                                            
      31 +    GC.gc(); t2=@elapsed begin O=O0; for _ in 1:steps; O=step(t
         +erms,O,dt,trunc,flow_end); O=trunc(O); end; println((end_len=le
         +ngth(O),)); end                                                
      32      println((N=N,M=M,steps=steps,stage=t1,endonly=t2))
      33  end
      34  bench(32,12,50)

  Ran 2 shell commands (ctrl+o to expand)

●

  Ran 1 shell command (ctrl+o to expand)

● Implemented a first generic TS-native Trotter prototype: TrotterTS.

  What it does

  Adds a new evolution method:

  TrotterTS(; order=2, substep=:taylor2)

  It works only for OperatorTS and does not expand into translated ordinary 
  gates.

  Instead, it splits the Hamiltonian into TS orbit terms:

  H = \sum_a H_a^{orb}

  and applies a Strang product over those orbit Liouvillians:

  e^{dt/2 L_{H_1^{orb}}}
  e^{dt/2 L_{H_2^{orb}}}
  ...
  e^{dt L_{H_m^{orb}}}
  ...
  e^{dt/2 L_{H_2^{orb}}}
  e^{dt/2 L_{H_1^{orb}}}

  Each orbit flow is currently approximated by a second-order Taylor action in
  the TS algebra:

  e^{dt L_a} O
  \approx
  O + dt L_a O + \frac{dt^2}{2} L_a^2 O

  with:

  L_a(O) = i[H_a^{orb}, O] / \hbar

  This is generic TS structure only — no MFIM commuting-layer/even-odd
  assumptions.

  Files changed

  - src/PauliStrings.jl
    - exported TrotterTS
  - src/evolution.jl
    - added TrotterTS method type
    - added orbit-term splitting
    - added TS orbit Liouvillian Taylor-2 subflow
    - added _evolve(::TrotterTS, H::OperatorTS, O::OperatorTS, ...)
  - test/evolution.jl
    - added small-system exact comparison for TrotterTS
    - checks output remains OperatorTS

  Also the earlier notes file is still present:

  benchmark/ts_trotter_slowdown_findings.md

  Correctness check

  Small N=6 exact comparison showed second-order scaling for TrotterTS:

  ┌────────┬─────────────────────────┐
  │   dt   │ relative error vs exact │
  ├────────┼─────────────────────────┤
  │    0.1 │                 3.02e-2 │
  ├────────┼─────────────────────────┤
  │   0.05 │                 7.50e-3 │
  ├────────┼─────────────────────────┤
  │  0.025 │                 1.87e-3 │
  ├────────┼─────────────────────────┤
  │ 0.0125 │                 4.68e-4 │
  └────────┴─────────────────────────┘

  Halving dt reduces error by ~4x, as expected for a second-order method.

  Benchmark smoke test

  MFIM, dt=0.1, no fout, same truncation style:

  ┌──────────────────────┬───────┬─────────────────┬───────────┐
  │         Case         │  RK4  │ current Trotter │ TrotterTS │
  ├──────────────────────┼───────┼─────────────────┼───────────┤
  │ N=32, M=12, steps=50 │ 1.45s │           1.43s │     1.10s │
  ├──────────────────────┼───────┼─────────────────┼───────────┤
  │ N=64, M=14, steps=20 │ 2.35s │           4.19s │     2.45s │
  └──────────────────────┴───────┴─────────────────┴───────────┘

  So the prototype already beats current expanded-gate Trotter, and is close to
  RK4 at N=64, faster at N=32. This confirms the direction is plausible, though
  it does not yet clearly outperform RK4 across the larger TS benchmark.

  Tests

  Full test suite passed:

  PauliStrings tests passed

  Command run:

  julia --project=/home/dhruva/Programs/unitaryHack/PauliStrings.jl -e 'using 
  Pkg; Pkg.test()'

  Important caveat

  This is a prototype TS-native second-order method, not an exact orbit
  exponential. The current subflow is Taylor-2, so accuracy/stability may be
  worse than an actual orbit exponential/Krylov subflow. The next meaningful
  step is likely to replace substep=:taylor2 with a Krylov/exponential-action
  orbit flow while preserving the same TS-orbit Strang structure.

✻ Baked for 6m 42s

❯ go ahead                                                                      

  Searched for 1 pattern (ctrl+o to expand)

● I’ll take the next step: add a Krylov/exponential-action substep for the same
  TS-orbit Strang structure, then compare it to the Taylor-2 substep and RK4.

  Searched for 1 pattern, read 1 file (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 7 lines, removed 4 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, substep=:taylor2)                     
      29 +    TrotterTS(; order=2, substep=:taylor2, krylovdim=4)        
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34 -Liouvillians. The current `substep=:taylor2` applies each orbit
         - flow with a                                                   
      35 -second-order Taylor approximation in the TS operator algebra.  
      34 +Liouvillians. `substep=:taylor2` applies each orbit flow with a
         + second-order Taylor                                           
      35 +approximation in the TS operator algebra. `substep=:krylov` app
         +lies a small Arnoldi                                           
      36 +exponential action with dimension `krylovdim`.                 
      37  """
      38  struct TrotterTS <: AbstractEvolutionMethod
      39      order::Int
      40      substep::Symbol
      41 +    krylovdim::Int                                             
      42  end
      41 -TrotterTS(; order::Integer=2, substep::Symbol=:taylor2) = Trott
         -erTS(Int(order), substep)                                      
      43 +TrotterTS(; order::Integer=2, substep::Symbol=:taylor2, krylovd
         +im::Integer=4) =                                               
      44 +    TrotterTS(Int(order), substep, Int(krylovdim))             
      45  
      46  """
      47      RK4()

● Update(src/evolution.jl)
  ⎿  Added 48 lines
      307      return truncation(O + dt * k1 + (dt^2 / 2) * k2)
      308  end
      309  
      310 +function _axpy_combination(basis, coeffs, template)           
      311 +    O = zero(template)                                        
      312 +    for j in eachindex(coeffs)                                
      313 +        iszero(coeffs[j]) || (O = O + basis[j] * coeffs[j])   
      314 +    end                                                       
      315 +    return O                                                  
      316 +end                                                           
      317 +                                                              
      318 +function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
          +Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
          +t, truncation)                                                
      319 +    krylovdim > 0 || throw(ArgumentError("krylovdim must be po
          +sitive"))                                                     
      320 +    β = sqrt(real(trace_product(O, O)))                       
      321 +    β == 0 && return O                                        
      322 +                                                              
      323 +    V = typeof(O)[]                                           
      324 +    push!(V, O / β)                                           
      325 +    Hk = zeros(ComplexF64, krylovdim + 1, krylovdim)          
      326 +    actual_dim = krylovdim                                    
      327 +                                                              
      328 +    for j in 1:krylovdim                                      
      329 +        w = truncation(_orbit_liouvillian(Ha, V[j], hbar))    
      330 +        for i in 1:j                                          
      331 +            Hk[i, j] = trace_product(V[i], w) / real(trace_pro
          +duct(V[i], V[i]))                                             
      332 +            w = truncation(w - V[i] * Hk[i, j])               
      333 +        end                                                   
      334 +        nw = sqrt(real(trace_product(w, w)))                  
      335 +        if nw < 1e-12                                         
      336 +            actual_dim = j                                    
      337 +            break                                             
      338 +        end                                                   
      339 +        Hk[j + 1, j] = nw                                     
      340 +        push!(V, w / nw)                                      
      341 +    end                                                       
      342 +                                                              
      343 +    F = exp(dt * Hk[1:actual_dim, 1:actual_dim])              
      344 +    coeffs = β * F[:, 1]                                      
      345 +    return truncation(_axpy_combination(V, coeffs, O))        
      346 +end                                                           
      347 +                                                              
      348 +function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
          +r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
          +runcation)                                                    
      349 +    if method.substep === :taylor2                            
      350 +        return _orbit_flow_taylor2(Ha, O, dt, hbar, truncation
          +)                                                             
      351 +    elseif method.substep === :krylov                         
      352 +        return _orbit_flow_krylov(Ha, O, dt, hbar, method.kryl
          +ovdim, truncation)                                            
      353 +    else                                                      
      354 +        throw(ArgumentError("unknown TrotterTS substep $(metho
          +d.substep); expected :taylor2 or :krylov"))                   
      355 +    end                                                       
      356 +end                                                           
      357 +                                                              
      358  function _orbit_terms(H::Operator{<:PauliStringTS})
      359      terms = typeof(H)[]
      360      for j in 1:length(H)

● Update(src/evolution.jl)
  ⎿  Added 4 lines, removed 5 lines
      366  
      367  function _trotterts_step(Hterms, O::Operator{<:PauliStringTS},
            dt::Real, method::TrotterTS, hbar::Real, truncation)
      368      method.order ∈ (1, 2) || throw(ArgumentError("order must b
           e 1 or 2, got $(method.order)"))
      369 -    method.substep === :taylor2 || throw(ArgumentError("only s
          -ubstep=:taylor2 is currently implemented"))                   
      369      isempty(Hterms) && return O
      370      if method.order == 1 || length(Hterms) == 1
      371          for Ha in Hterms
      373 -            O = _orbit_flow_taylor2(Ha, O, dt, hbar, truncatio
          -n)                                                            
      372 +            O = _orbit_flow(Ha, O, dt, method, hbar, truncatio
          +n)                                                            
      373          end
      374      else
      375          L = length(Hterms)
      376          for j in 1:(L - 1)
      378 -            O = _orbit_flow_taylor2(Hterms[j], O, dt / 2, hbar
          -, truncation)                                                 
      377 +            O = _orbit_flow(Hterms[j], O, dt / 2, method, hbar
          +, truncation)                                                 
      378          end
      380 -        O = _orbit_flow_taylor2(Hterms[L], O, dt, hbar, trunca
          -tion)                                                         
      379 +        O = _orbit_flow(Hterms[L], O, dt, method, hbar, trunca
          +tion)                                                         
      380          for j in (L - 1):-1:1
      382 -            O = _orbit_flow_taylor2(Hterms[j], O, dt / 2, hbar
          -, truncation)                                                 
      381 +            O = _orbit_flow(Hterms[j], O, dt / 2, method, hbar
          +, truncation)                                                 
      382          end
      383      end
      384      return O

  Ran 2 shell commands (ctrl+o to expand)
● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 3 lines
      317  
      318  function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
           t, truncation)
      319      krylovdim > 0 || throw(ArgumentError("krylovdim must be po
           sitive"))
      320 -    β = sqrt(real(trace_product(O, O)))                       
      320 +    β = sqrt(real(trace_product(O', O)))                      
      321      β == 0 && return O
      322  
      323      V = typeof(O)[]
     ...
      328      for j in 1:krylovdim
      329          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      330          for i in 1:j
      331 -            Hk[i, j] = trace_product(V[i], w) / real(trace_pro
          -duct(V[i], V[i]))                                             
      331 +            Hk[i, j] = trace_product(V[i]', w) / trace_product
          +(V[i]', V[i])                                                 
      332              w = truncation(w - V[i] * Hk[i, j])
      333          end
      334 -        nw = sqrt(real(trace_product(w, w)))                  
      334 +        nw = sqrt(real(trace_product(w', w)))                 
      335          if nw < 1e-12
      336              actual_dim = j
      337              break

  Read 1 file, ran 3 shell commands (ctrl+o to expand)

● Update(test/evolution.jl)
  ⎿  Added 2 lines
      33      res4 = evolve(resum(H), resum(O0), times; method = Trotter(
          ), fout=fout_trotter, truncation=truncation).history
      34      res5 = evolve(H, O0, times; method = Trotter(), fout=fout, 
          truncation=truncation).history
      35      res6 = evolve(H, O0, times; method = TrotterTS(), fout=fout
          , truncation=truncation).history
      36 +    res7 = evolve(H, O0, times; method = TrotterTS(substep=:kry
         +lov, krylovdim=4), fout=fout, truncation=truncation).history   
      37  
      38      @test norm(res2 .- res1)/norm(res1) < 1e-7
      39      @test norm(res3 .- res1)/norm(res1) < 1e-7
      40      @test norm(res4 .- res1)/norm(res1) < 1e-5
      41      @test norm(res5 .- res1)/norm(res1) < 1e-4
      42      @test norm(res6 .- res1)/norm(res1) < 1e-3
      43 +    @test norm(res7 .- res1)/norm(res1) < 1e-3                 
      44      @test evolve(H, O0, [0, dt]; method = TrotterTS(), fout=not
          hing, truncation=truncation).final isa typeof(O0)
      45  
      46      dissipation(O, dt) = add_noise(O, 0.1*dt)

  Ran 2 shell commands (ctrl+o to expand)

● Added the next substep option: TrotterTS(substep=:krylov, krylovdim=...).

  What changed

  TrotterTS now supports:

  TrotterTS(; order=2, substep=:taylor2, krylovdim=4)
  TrotterTS(; order=2, substep=:krylov, krylovdim=4)

  The structure is the same TS-orbit Strang method:

  H = \sum_a H_a^{orb}

  but for each orbit flow:

  e^{\tau L_{H_a^{orb}}} O

  we can now use a small Arnoldi/Krylov exponential action instead of Taylor-2.

  I used the TS Hilbert-Schmidt inner product via:

  trace_product(V[i]', w)

  so the Krylov basis handles complex coefficients correctly.

  Accuracy check

  Small N=6, exact dense comparison, TFIM, relative trajectory error:

  TrotterTS(substep=:taylor2)

  ┌────────┬─────────┐
  │   dt   │  error  │
  ├────────┼─────────┤
  │    0.1 │ 3.02e-2 │
  ├────────┼─────────┤
  │   0.05 │ 7.50e-3 │
  ├────────┼─────────┤
  │  0.025 │ 1.87e-3 │
  ├────────┼─────────┤
  │ 0.0125 │ 4.68e-4 │
  └────────┴─────────┘

  Second-order scaling, but fairly large constant.

  TrotterTS(substep=:krylov, krylovdim=4)

  ┌───────┬─────────┐
  │  dt   │  error  │
  ├───────┼─────────┤
  │   0.1 │ 1.44e-3 │
  ├───────┼─────────┤
  │  0.05 │ 3.71e-4 │
  ├───────┼─────────┤
  │ 0.025 │ 9.46e-5 │
  └───────┴─────────┘

  Also second-order scaling, but much better accuracy constant. krylovdim=6 did
  not improve these numbers on this test, so the remaining error is likely the
  outer orbit-Strang splitting, not the orbit exponential approximation.

  krylovdim=2 was poor:

  ┌───────┬─────────┐
  │  dt   │  error  │
  ├───────┼─────────┤
  │   0.1 │ 1.80e-1 │
  ├───────┼─────────┤
  │  0.05 │ 9.22e-2 │
  ├───────┼─────────┤
  │ 0.025 │ 4.65e-2 │
  └───────┴─────────┘

  Performance smoke benchmark

  MFIM, dt=0.1, same truncation, no fout.

  ┌──────────┬───────┬─────────┬─────────────────────┬──────────────────────┐
  │          │       │ current │                     │                      │
  │   Case   │  RK4  │         │ TrotterTS(:taylor2) │ TrotterTS(:krylov,4) │
  │          │       │ Trotter │                     │                      │
  ├──────────┼───────┼─────────┼─────────────────────┼──────────────────────┤
  │ N=32,    │       │         │                     │                      │
  │ M=12,    │ 1.42s │   1.59s │               1.09s │                5.56s │
  │ steps=50 │       │         │                     │                      │
  ├──────────┼───────┼─────────┼─────────────────────┼──────────────────────┤
  │ N=64,    │       │         │                     │                      │
  │ M=14,    │ 2.33s │   4.33s │               2.44s │                6.30s │
  │ steps=20 │       │         │                     │                      │
  └──────────┴───────┴─────────┴─────────────────────┴──────────────────────┘

  So:

  - :krylov improves accuracy a lot, but is too slow in the current
  straightforward implementation.
  - :taylor2 is faster and TS-native, but much less accurate.
  - The TS-orbit decomposition direction is still validated, but the subflow
  implementation matters.

  Tests

  Full package tests pass:

  PauliStrings tests passed

  Current state

  Changed files:

  src/PauliStrings.jl
  src/evolution.jl
  test/evolution.jl
  benchmark/ts_trotter_slowdown_findings.md

  The important conclusion from this step: exact-ish orbit exponentials via 
  naïve Krylov are accurate but too expensive; Taylor-2 is fast but lower 
  accuracy. The next practical avenue is probably a better low-cost orbit
  subflow or optimizing the TS orbit commutator kernel, not increasing Krylov
  dimension.

✻ Worked for 5m 26s

❯ before you continue, this must be an explored algorithm/problem in            
  literature. Look at past research, and then proceed.                          

● Agreed — before pushing the algorithm further, we should anchor it in known
  work. I’ll pause implementation and do a literature pass specifically for
  product formulas / time evolution exploiting translation invariance or 
  symmetry-reduced operator spaces, then map what applies to this Pauli-string
  TS setting.

● Fetch(https://arxiv.org/abs/math-ph/0506007)
  ⎿  Received 45.3KB (200 OK)

● Fetch(https://arxiv.org/abs/1810.07006)
  ⎿  Received 46.6KB (200 OK)

● Fetch(https://arxiv.org/abs/1008.3477)
  ⎿  Received 46.6KB (200 OK)

● Fetch(https://arxiv.org/abs/1007.4145)
  ⎿  Received 47.5KB (200 OK)

● Web Search("translation invariant Hamiltonian simulation Trotter product
            formula symmetry reduced subspace")
  ⎿  Did 0 searches in 673ms

● Agreed. I paused implementation and did a literature pass. The short version
  is:

  ▎ There is strong literature support for (1) product formulas, (2) 
  ▎ Krylov/exponential action, and (3) translation-invariant time evolution in 
  ▎ compressed manifolds such as iTEBD/uMPS/TDVP. But I did not find a directly 
  ▎ standard “translation-symmetric Pauli-string orbit Trotter primitive” 
  ▎ matching this repo’s OperatorTS representation. Our design is closest in 
  ▎ spirit to splitting the dynamics over symmetry-preserving orbit 
  ▎ Liouvillians, analogous to tangent-space / infinite-system tensor-network 
  ▎ methods, but implemented in a Pauli-string TS operator basis.

  Relevant literature buckets

  1. Product formulas / Trotter-Suzuki splitting

  This is the general foundation.

  - Hatano & Suzuki, “Finding Exponential Product Formulas of Higher Orders”
  reviews Trotter/Suzuki product formulas and higher-order decompositions. It
  explicitly frames product formulas as preserving important dynamical
  symmetries.
  Source: arXiv:math-ph/0506007
  - Childs et al., “A Theory of Trotter Error” gives modern error analysis for
  product formulas, especially exploiting locality and commutativity structure.
  It is relevant for understanding why product formulas can be efficient for
  local Hamiltonians, but it does not appear to focus on translation-invariant
  compressed representations.
  Source: arXiv:1912.08854
  - Haah, Hastings, Kothari & Low, “Quantum algorithm for simulating real time
  evolution of lattice Hamiltonians” gives near-optimal lattice Hamiltonian
  simulation algorithms based on locality/Lieb-Robinson ideas and products of
  small unitaries. Again, locality is central; translation-invariant compression
  is not the main axis.
  Source: arXiv:1801.03922

  Implication for us:
  Standard Trotter theory justifies Strang splitting, but current TS Trotter’s
  problem is the choice of summands. Splitting into individual translated Pauli
  terms gives O(N) factors. Splitting into translation-orbit Hamiltonians is the
  symmetry-preserving alternative.

  ---
  2. Translation-invariant time evolution in compressed representations

  This is the closest conceptual literature.

  - Vidal’s infinite TEBD paper explicitly exploits translation invariance to
  simulate infinite 1D quantum lattice systems and time evolution. It works in
  an infinite tensor-network representation rather than a finite expanded
  lattice.
  Source: arXiv:cond-mat/0605597
  - McCulloch’s infinite-size DMRG work targets fixed-point
  translation-invariant MPS wavefunctions and is related to infinite TEBD.
  Source: arXiv:0804.2509
  - Haegeman et al., “Unifying time evolution and optimization with matrix
  product states” formulates time evolution using TDVP/tangent-space projection
  rather than ordinary gate decomposition.
  Source: arXiv:1408.5056
  - Vanderstraeten, Haegeman & Verstraete’s uniform MPS tangent-space lecture
  notes cover tangent spaces and real-time evolution for uniform MPS, i.e.
  translation-invariant thermodynamic-limit representations.
  Source: arXiv:1810.07006

  Implication for us:
  The design principle is well supported:

  do not expand the translation-invariant system into all sites;
  evolve directly inside the translation-invariant/compressed manifold.

  But those works use MPS/uMPS tangent-space geometry, not Pauli-string
  OperatorTS.

  Our TrotterTS is more like a symmetry-reduced operator-basis analogue:

  uMPS/TDVP:      project dynamics to translation-invariant MPS manifold
  OperatorTS:     keep dynamics in translation-invariant Pauli-string orbit
  space

  ---
  3. Krylov/exponential-action methods

  The Krylov substep we added is standard numerically.

  - Krylov/Arnoldi approximations to exp(tA)v are a classic way to apply matrix
  exponentials without forming the matrix. I attempted to fetch Saad’s 1992 SIAM
  paper via DOI, but the fetch returned a different SIAM page, so I do not have
  source-supported metadata from the tool. Still, the method itself is standard
  and our usage is conventional: matrix-free action of an orbit Liouvillian.

  Implication for us:
  TrotterTS(substep=:krylov) is not a weird algorithmic invention; it is
  standard exponential-action machinery applied to a TS orbit Liouvillian.
  However, our naive implementation is too slow in current benchmarks.

  ---
  4. Pauli propagation / Pauli-string evolution

  - Rudolph et al., “Pauli Propagation: A Computational Framework for Simulating
  Quantum Systems” is directly relevant to Pauli-string/operator evolution and
  bit-level Pauli simulation. The arXiv page says it covers Pauli propagation as
  a classical simulation framework for circuits and dynamics, with a Julia
  package. The page did not expose details on translation symmetry or product
  formulas in the fetch, but it is relevant background for Pauli-string
  propagation.
  Source: arXiv:2505.21606

  Implication for us:
  The existing non-TS trotter_step! is in this Pauli propagation spirit: single
  Pauli rotations map one Pauli string to at most two branches. The missing
  piece is a translation-invariant/orbit-level analogue.

  ---
  What this means for our design

  The strongest literature-backed route is:

  Strang splitting over translation-orbit Hamiltonians
  +
  matrix-free action of each orbit Liouvillian inside OperatorTS

  Formally:

  H = \sum_a H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_r T_r P_a T_r^\dagger

  and:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2}L_{H_1^{orb}}}
  \cdots
  e^{\Delta t L_{H_m^{orb}}}
  \cdots
  e^{\frac{\Delta t}{2}L_{H_1^{orb}}}.

  This is exactly analogous to respecting translation symmetry at the
  algorithmic level, as iTEBD/uMPS/TDVP do, but in the Pauli-string operator
  basis.

  Important caution

  The Taylor-2 orbit substep is not really a standard “Trotter subflow”; it is a
  cheap second-order approximation to each orbit exponential. It is useful as a
  prototype, but we should describe it carefully:

  TrotterTS(:taylor2) = orbit-Strang with Taylor-2 approximate orbit flows

  The Krylov orbit substep is more literature-standard:

  TrotterTS(:krylov) = orbit-Strang with matrix-free exponential action

  but our current naive version is too slow.

  Where to proceed

  Given the literature, I would not continue blindly with Taylor-2 as the
  “final” method. I’d treat it as a fast prototype showing that orbit-level TS
  splitting has potential.

  The next research-aligned step should be:

  1. Keep orbit-level Strang as the core method.
  2. Improve the matrix-free orbit Liouvillian application:
    - the current TS commutator scans all shifts;
    - for local Pauli strings, many shifts cannot contribute;
    - support-aware relative-shift enumeration is the right optimization.
  3. Re-test :krylov after the orbit Liouvillian is faster.
  4. Compare at matched accuracy against RK4.

  That is the most defensible literature-backed path: TS orbit splitting + 
  matrix-free exponential action, optimized using locality of Pauli support but
  not MFIM-specific commuting structure.

✻ Churned for 3m 25s

❯ Ok, continue. Please remember that Pauli Propogation is exactly the package   
  that depends on PauliStrings.jl which was written by the same authors. So, it 
  is doubtful that you will get the solution there, but you will find good      
  context. When doing performance tests, also check if the slowdowns are due to 
  memory allocations rather than algorithmic slowdowns. However, do not get     
  rabbit holed trying to do this. again, test that the reduction in the number  
  of terms is actually happening. Do not blindly trust the final time number,   
  instead actually check.                                                       

  Recalled 1 memory, wrote 2 memories (ctrl+o to expand)

● Got it. I’ll continue, but with two guardrails:

  1. Treat PauliPropagation as context only, not a likely source of the
  solution.
  2. For every performance claim, also inspect allocations and intermediate term
  counts, so we verify the TS-orbit method is actually reducing work rather
  than just moving costs around.

● Write(~/.claude/jobs/513a65b3/tmp/trotterts_instrument.jl)
  ⎿  Wrote 152 lines to
     ../../../.claude/jobs/513a65b3/tmp/trotterts_instrument.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N, h)
       5     H = Operator(N)
       6     H += -h, "X", 1
       7     H += -h/2, "Z", 2
       8     H += "Z", 1, "Z", 2
       9     return OperatorTS{(N,)}(H)
      10 end
     … +161 lines (ctrl+o to expand)
  Ran 1 shell command (ctrl+o to expand)

● I continued with instrumentation rather than more algorithm changes.

  I also saved your PauliPropagation note as memory: PauliPropagation.jl is
  useful context, but not likely where the solution lives because it depends on
  PauliStrings.jl and has the same author lineage.

  Instrumentation result

  I checked, for one representative warmed step, not just final runtime:

  - number of expanded gates vs orbit factors,
  - Liouvillian/commutator calls,
  - truncation calls,
  - max intermediate term count,
  - allocations,
  - time spent in commutators/truncation.

  Script used:

  /home/dhruva/.claude/jobs/513a65b3/tmp/trotterts_instrument.jl

  N=32, M=12

  lenHts = 3
  lenHr = 96
  expanded_gates = 191
  orbit_factors = 5
  start_len = 4096

  ┌────────────────────┬────────┬────────────┬─────────────┬──────┬─────┐
  │                    │        │ Allocation │    Main     │ Trun │ Max │
  │       Method       │  Time  │     s      │ operation   │ c ca │     │
  │                    │        │            │    count    │ lls  │ len │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ RK4                │ 0.0189 │    18.1 MB │       4 TS  │    5 │ 763 │
  │                    │      s │            │ commutators │      │   0 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ Expanded Trotter   │ 0.0432 │    57.2 MB │   191 gate  │  192 │ 614 │
  │                    │      s │            │     updates │      │   7 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ TrotterTS(:taylor2 │ 0.0216 │    28.6 MB │ 10 orbit Li │   16 │ 642 │
  │ )                  │      s │            │  ouvillians │      │   9 │
  └────────────────────┴────────┴────────────┴─────────────┴──────┴─────┘

  N=64, M=14

  lenHts = 3
  lenHr = 192
  expanded_gates = 383
  orbit_factors = 5
  start_len = 16384

  ┌────────────────────┬────────┬────────────┬─────────────┬──────┬──────┐
  │                    │        │ Allocation │    Main     │ Trun │ Max  │
  │       Method       │  Time  │     s      │ operation   │ c ca │ len  │
  │                    │        │            │    count    │ lls  │      │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼──────┤
  │ RK4                │ 0.1358 │    86.2 MB │       4 TS  │    5 │ 3003 │
  │                    │      s │            │ commutators │      │    8 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼──────┤
  │ Expanded Trotter   │ 0.2754 │   470.7 MB │   383 gate  │  384 │ 2466 │
  │                    │      s │            │     updates │      │    9 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼──────┤
  │ TrotterTS(:taylor2 │ 0.1290 │   138.9 MB │ 10 orbit Li │   16 │ 2637 │
  │ )                  │      s │            │  ouvillians │      │    3 │
  └────────────────────┴────────┴────────────┴─────────────┴──────┴──────┘

  What this verifies

  The structural reduction is real:

  N=64 expanded Trotter: 383 ordinary gate updates
  N=64 TrotterTS:          5 orbit factors

  For Taylor-2 each orbit factor costs 2 Liouvillian applications, so:

  5 orbit factors -> 10 TS Liouvillian calls

  That is happening as expected.

  The truncation reduction is also real:

  N=64 expanded Trotter: 384 truncation calls
  N=64 TrotterTS:         16 truncation calls

  Allocations also confirm this is not just moving the same work around:

  N=64 expanded Trotter: ~471 MB
  N=64 TrotterTS:        ~139 MB
  RK4:                   ~86 MB

  So TrotterTS(:taylor2) substantially reduces work and memory relative to
  expanded Trotter, but it still allocates more than RK4.

  Current diagnosis

  The bottleneck after moving to TS-orbit splitting is no longer the expanded
  gate count. It is now the cost of each TS orbit Liouvillian application:

  TrotterTS(:taylor2) N=64:
    total time:     0.129s
    commutator time 0.103s

  So ~80% of the time is inside TS commutators.

  That means the next real target is not more Trotter orchestration; it is the
  TS orbit commutator kernel:

  commutator(Ha_orbit, Ots)

  At present, this still uses the generic TS binary kernel, which loops over
  translations. The algorithmic reduction from 383 ordinary gates to 10 orbit
  Liouvillian calls is real, but each orbit Liouvillian is still fairly
  expensive.

  Important comparison

  TrotterTS(:krylov, krylovdim=4) improved accuracy a lot but was slow:

  ┌─────────────────────┬───────┬─────────────────────┬──────────────────────┐
  │        Case         │  RK4  │ TrotterTS(:taylor2) │ TrotterTS(:krylov,4) │
  ├─────────────────────┼───────┼─────────────────────┼──────────────────────┤
  │ N=32, M=12,         │ 1.42s │               1.09s │                5.56s │
  │ steps=50            │       │                     │                      │
  ├─────────────────────┼───────┼─────────────────────┼──────────────────────┤
  │ N=64, M=14,         │ 2.33s │               2.44s │                6.30s │
  │ steps=20            │       │                     │                      │
  └─────────────────────┴───────┴─────────────────────┴──────────────────────┘

  So we should not continue optimizing Krylov first. It is accurate but too
  allocation/commutator-heavy in the current implementation.

  Next step I recommend

  Keep TrotterTS(:taylor2) as the fast proof-of-concept and optimize/measure the
  TS orbit Liouvillian:

  commutator(single_orbit_term, Ots)

  Specifically, check whether the current kernel spends time scanning relative
  shifts that cannot contribute. If yes, implement a support-aware
  relative-shift enumeration for Pauli strings, then re-run the same
  instrumentation:

  - operation counts,
  - max term counts,
  - allocations,
  - commutator time,
  - total time.

  This is still generic TS/Pauli-string structure, not MFIM-specific commuting
  structure.

✻ Crunched for 2m 39s

❯ Continue                                                                      

● I’ll continue by targeting the TS orbit Liouvillian bottleneck: add a
  support-aware relative-shift path for TS commutators/products, then verify it
  reduces shift work, allocations, and time rather than only improving final
  timings.

  Read 3 files, ran 2 shell commands (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Added 34 lines
      123  end
      124  @inline all_shifts(::Type{<:PauliStringTS{Ls,Ps}}) where {Ls,P
           s} = all_shifts(Ls, Ps)
      125  
      126 +function _site_coords(i::Integer, Ls::Tuple)                  
      127 +    r = i - 1                                                 
      128 +    return ntuple(length(Ls)) do k                            
      129 +        L = Ls[k]                                             
      130 +        c = mod(r, L) + 1                                     
      131 +        r = div(r, L)                                         
      132 +        c                                                     
      133 +    end                                                       
      134 +end                                                           
      135 +                                                              
      136 +function _overlap_shifts(p1::PauliString, p2::PauliString, Ls:
          +:Tuple, Ps::Tuple)                                            
      137 +    s1 = support(p1)                                          
      138 +    s2 = support(p2)                                          
      139 +    (isempty(s1) || isempty(s2)) && return typeof(first(all_sh
          +ifts(Ls, Ps)))[]                                              
      140 +                                                              
      141 +    shifts = Set{typeof(first(all_shifts(Ls, Ps)))}()         
      142 +    coords1 = _site_coords.(s1, Ref(Ls))                      
      143 +    coords2 = _site_coords.(s2, Ref(Ls))                      
      144 +    for c1 in coords1, c2 in coords2                          
      145 +        valid = true                                          
      146 +        shift_tuple = ntuple(length(Ls)) do k                 
      147 +            if Ps[k]                                          
      148 +                d = mod(c1[k] - c2[k], Ls[k])                 
      149 +                iszero(d) ? Ls[k] : d                         
      150 +            else                                              
      151 +                valid &= c1[k] == c2[k]                       
      152 +                1                                             
      153 +            end                                               
      154 +        end                                                   
      155 +        valid && push!(shifts, shift_tuple)                   
      156 +    end                                                       
      157 +    return collect(shifts)                                    
      158 +end                                                           
      159 +                                                              
      160  """
      161      shift(p::PauliString, Ls::Tuple, shifts::Tuple)
      162      shift(p::PauliString, Ls::Tuple, Ps::Tuple, shifts::Tuple)

● Update(src/operations_strings.jl)
  ⎿  Added 1 line, removed 1 line
      141      @inbounds for (p1, c1) in zip(p1s, c1s)
      142          rep1 = representative(p1)
      143          rep2 = representative(p2)
      144 -        for s in all_shifts(paulistringtype(A))               
      144 +        for s in _overlap_shifts(rep1, rep2, Ls, Ps)          
      145              p, k = op(rep1, shift(rep2, Ls, Ps, s))
      146              c = c1 * c2 * k
      147              if (k != 0) && (abs(c) > epsilon) && pauli_weight(
           p) < maxlength

● Update(src/operations_strings.jl)
  ⎿  Added 1 line, removed 1 line
      184      # core kernel logic
      185      rep1 = representative(p1)
      186      rep2 = representative(p2)
      187 -    for s in all_shifts(paulistringtype(A))                   
      187 +    for s in _overlap_shifts(rep1, rep2, Ls, Ps)              
      188          p, k = op(rep1, shift(rep2, Ls, Ps, s))
      189          c = c1 * c2 * k
      190          if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) <
            maxlength

  Read 1 file (ctrl+o to expand)

● Update(src/operations_strings.jl)
  ⎿  Added 1 line, removed 1 line
      141      @inbounds for (p1, c1) in zip(p1s, c1s)
      142          rep1 = representative(p1)
      143          rep2 = representative(p2)
      144 -        for s in _overlap_shifts(rep1, rep2, Ls, Ps)          
      144 +        for s in all_shifts(paulistringtype(A))               
      145              p, k = op(rep1, shift(rep2, Ls, Ps, s))
      146              c = c1 * c2 * k
      147              if (k != 0) && (abs(c) > epsilon) && pauli_weight(
           p) < maxlength

● Update(src/translation_symmetry.jl)
  ⎿  Added 2 lines, removed 1 line
      309          rep1 = representative(p1)
      310          for (p2, c2) in zip(p2s, c2s)
      311              rep2 = representative(p2)
      312 -            for s in all_shifts(paulistringtype(A))           
      312 +            shifts = op === commutator ? _overlap_shifts(rep1,
          + rep2, Ls, Ps) : all_shifts(paulistringtype(A))               
      313 +            for s in shifts                                   
      314                  p, k = op(rep1, shift(rep2, Ls, Ps, s))
      315                  c = c1 * c2 * k
      316                  if (k != 0) && pauli_weight(p) < maxlength

  Read 1 file, ran 3 shell commands (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Added 30 lines, removed 23 lines
      123  end
      124  @inline all_shifts(::Type{<:PauliStringTS{Ls,Ps}}) where {Ls,P
           s} = all_shifts(Ls, Ps)
      125  
      126 -function _site_coords(i::Integer, Ls::Tuple)                  
      126 +@inline function _site_coord(i::Integer, Ls::Tuple, k::Integer
          +)                                                             
      127      r = i - 1
      128 -    return ntuple(length(Ls)) do k                            
      129 -        L = Ls[k]                                             
      130 -        c = mod(r, L) + 1                                     
      131 -        r = div(r, L)                                         
      132 -        c                                                     
      128 +    for j in 1:(k - 1)                                        
      129 +        r = div(r, Ls[j])                                     
      130      end
      131 +    return mod(r, Ls[k]) + 1                                  
      132  end
      133  
      134  function _overlap_shifts(p1::PauliString, p2::PauliString, Ls:
           :Tuple, Ps::Tuple)
      137 -    s1 = support(p1)                                          
      138 -    s2 = support(p2)                                          
      139 -    (isempty(s1) || isempty(s2)) && return typeof(first(all_sh
          -ifts(Ls, Ps)))[]                                              
      135 +    u1 = p1.v | p1.w                                          
      136 +    u2 = p2.v | p2.w                                          
      137 +    (iszero(u1) || iszero(u2)) && return typeof(first(all_shif
          +ts(Ls, Ps)))[]                                                
      138  
      141 -    shifts = Set{typeof(first(all_shifts(Ls, Ps)))}()         
      142 -    coords1 = _site_coords.(s1, Ref(Ls))                      
      143 -    coords2 = _site_coords.(s2, Ref(Ls))                      
      144 -    for c1 in coords1, c2 in coords2                          
      145 -        valid = true                                          
      146 -        shift_tuple = ntuple(length(Ls)) do k                 
      147 -            if Ps[k]                                          
      148 -                d = mod(c1[k] - c2[k], Ls[k])                 
      149 -                iszero(d) ? Ls[k] : d                         
      150 -            else                                              
      151 -                valid &= c1[k] == c2[k]                       
      152 -                1                                             
      139 +    shifts = typeof(first(all_shifts(Ls, Ps)))[]              
      140 +    seen = Set{typeof(first(all_shifts(Ls, Ps)))}()           
      141 +    N = qubitlength(p1)                                       
      142 +    for i1 in 1:N                                             
      143 +        ((u1 >> (i1 - 1)) & 1 == 0) && continue               
      144 +        for i2 in 1:N                                         
      145 +            ((u2 >> (i2 - 1)) & 1 == 0) && continue           
      146 +            valid = true                                      
      147 +            shift_tuple = ntuple(length(Ls)) do k             
      148 +                c1 = _site_coord(i1, Ls, k)                   
      149 +                c2 = _site_coord(i2, Ls, k)                   
      150 +                if Ps[k]                                      
      151 +                    d = mod(c1 - c2, Ls[k])                   
      152 +                    iszero(d) ? Ls[k] : d                     
      153 +                else                                          
      154 +                    valid &= c1 == c2                         
      155 +                    1                                         
      156 +                end                                           
      157              end
      158 +            if valid && !(shift_tuple in seen)                
      159 +                push!(seen, shift_tuple)                      
      160 +                push!(shifts, shift_tuple)                    
      161 +            end                                               
      162          end
      155 -        valid && push!(shifts, shift_tuple)                   
      163      end
      157 -    return collect(shifts)                                    
      164 +    return shifts                                             
      165  end
      166  
      167  """

  Ran 2 shell commands (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Removed 41 lines
      123  end
      124  @inline all_shifts(::Type{<:PauliStringTS{Ls,Ps}}) where {Ls,P
           s} = all_shifts(Ls, Ps)
      125  
      126 -@inline function _site_coord(i::Integer, Ls::Tuple, k::Integer
          -)                                                             
      127 -    r = i - 1                                                 
      128 -    for j in 1:(k - 1)                                        
      129 -        r = div(r, Ls[j])                                     
      130 -    end                                                       
      131 -    return mod(r, Ls[k]) + 1                                  
      132 -end                                                           
      133 -                                                              
      134 -function _overlap_shifts(p1::PauliString, p2::PauliString, Ls:
          -:Tuple, Ps::Tuple)                                            
      135 -    u1 = p1.v | p1.w                                          
      136 -    u2 = p2.v | p2.w                                          
      137 -    (iszero(u1) || iszero(u2)) && return typeof(first(all_shif
          -ts(Ls, Ps)))[]                                                
      138 -                                                              
      139 -    shifts = typeof(first(all_shifts(Ls, Ps)))[]              
      140 -    seen = Set{typeof(first(all_shifts(Ls, Ps)))}()           
      141 -    N = qubitlength(p1)                                       
      142 -    for i1 in 1:N                                             
      143 -        ((u1 >> (i1 - 1)) & 1 == 0) && continue               
      144 -        for i2 in 1:N                                         
      145 -            ((u2 >> (i2 - 1)) & 1 == 0) && continue           
      146 -            valid = true                                      
      147 -            shift_tuple = ntuple(length(Ls)) do k             
      148 -                c1 = _site_coord(i1, Ls, k)                   
      149 -                c2 = _site_coord(i2, Ls, k)                   
      150 -                if Ps[k]                                      
      151 -                    d = mod(c1 - c2, Ls[k])                   
      152 -                    iszero(d) ? Ls[k] : d                     
      153 -                else                                          
      154 -                    valid &= c1 == c2                         
      155 -                    1                                         
      156 -                end                                           
      157 -            end                                               
      158 -            if valid && !(shift_tuple in seen)                
      159 -                push!(seen, shift_tuple)                      
      160 -                push!(shifts, shift_tuple)                    
      161 -            end                                               
      162 -        end                                                   
      163 -    end                                                       
      164 -    return shifts                                             
      165 -end                                                           
      166 -                                                              
      126  """
      127      shift(p::PauliString, Ls::Tuple, shifts::Tuple)
      128      shift(p::PauliString, Ls::Tuple, Ps::Tuple, shifts::Tuple)

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 2 lines
      275          rep1 = representative(p1)
      276          for (p2, c2) in zip(p2s, c2s)
      277              rep2 = representative(p2)
      278 -            shifts = op === commutator ? _overlap_shifts(rep1,
          - rep2, Ls, Ps) : all_shifts(paulistringtype(A))               
      279 -            for s in shifts                                   
      278 +            for s in all_shifts(paulistringtype(A))           
      279                  p, k = op(rep1, shift(rep2, Ls, Ps, s))
      280                  c = c1 * c2 * k
      281                  if (k != 0) && pauli_weight(p) < maxlength

  Ran 2 shell commands (ctrl+o to expand)

● Update(src/operations_strings.jl)
  ⎿  Added 1 line, removed 1 line
      184      # core kernel logic
      185      rep1 = representative(p1)
      186      rep2 = representative(p2)
      187 -    for s in _overlap_shifts(rep1, rep2, Ls, Ps)              
      187 +    for s in all_shifts(paulistringtype(A))                   
      188          p, k = op(rep1, shift(rep2, Ls, Ps, s))
      189          c = c1 * c2 * k
      190          if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) <
            maxlength
  Ran 1 shell command (ctrl+o to expand)

● I continued with the support-aware commutator idea, but the first attempt was
  a regression, so I reverted it.

  What I tried

  I added an experimental _overlap_shifts helper to avoid scanning all
  translations in the TS commutator. The idea was:

  only shifts where support(rep1) overlaps support(shift(rep2))

  Then use this only for commutator(A::OperatorTS, B::OperatorTS), since
  zero-overlap shifts always commute and contribute nothing.

  Correctness spot check passed:

  norm(resum(commutator(Hts, Ots)) - commutator(resum(Hts), resum(Ots))) == 0

  for small N.

  Why I reverted it

  Performance got much worse.

  Before support-aware shift filtering, instrumentation showed roughly:

  N=64 TrotterTS-Taylor2:
    time  ≈ 0.129s
    alloc ≈ 139 MB
    comm_time ≈ 0.103s

  After the first support-aware version:

  N=64 TrotterTS-Taylor2:
    time  ≈ 0.676s
    alloc ≈ 665 MB
    comm_time ≈ 0.638s

  Even an optimized version was still worse:

  N=64 TrotterTS-Taylor2:
    time  ≈ 0.254s
    alloc ≈ 457 MB
    comm_time ≈ 0.182s

  So although the number of candidate shifts was reduced, the overhead of
  computing support overlaps, allocating shift vectors/sets, and doing per-pair
  bookkeeping dominated.

  I reverted the overlap-shift changes. test/evolution.jl passes again.

  Current state

  Kept:

  - TrotterTS(:taylor2)
  - TrotterTS(:krylov)
  - tests for both
  - documentation/notes file

  Reverted:

  - support-aware _overlap_shifts commutator optimization

  Lesson

  The “support-aware shifts” idea is still plausible, but it cannot be
  implemented with per-pair dynamic allocation. A viable version would need to
  be allocation-free or cached.

  A better next attempt would be one of:

  1. Cache active shifts per Hamiltonian orbit term and Pauli support pattern,
  not recompute per pair.
  2. Specialize for 1D contiguous/local supports without allocating
  sets/vectors.
  3. Avoid support lists entirely, using bit masks and shifted masks to check
  overlap cheaply.
  4. Instrument shift counts without changing execution, to estimate theoretical
  savings before another implementation attempt.

  Given your warning about not rabbit-holing on allocations: I’d do option 4
  next — measure how many full shifts are actually nonzero contributors in
  current TS commutator, without changing the algorithm, then decide whether a
  low-allocation cached implementation is worth it.

✻ Worked for 7m 5s

❯ continue                                                                      

● I’ll continue by measuring the actual shift-level sparsity in the TS orbit
  commutators, and also check whether the current kernel is doing unnecessary
  expensive shifts of the large O string rather than the small Hamiltonian orbit
  representative.

● Write(~/.claude/jobs/513a65b3/tmp/shift_sparsity_measure.jl)
  ⎿  Wrote 88 lines to
     ../../../.claude/jobs/513a65b3/tmp/shift_sparsity_measure.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N, h)
       5     H = Operator(N)
       6     H += -h, "X", 1
       7     H += -h/2, "Z", 2
       8     H += "Z", 1, "Z", 2
       9     return OperatorTS{(N,)}(H)
      10 end
     … +84 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/shift_sparsity_measure.jl)
  ⎿  Added 1 line, removed 1 line
      39          rep2 = representative(p2)
      40          w = pauli_weight(rep2)
      41          totw, nzw = get(len_by_weight, w, (0,0))
      42 -        for s in all_shifts(paulistringtype(O))                
      42 +        for s in PauliStrings.all_shifts(paulistringtype(O))   
      43              total += 1
      44              totw += 1
      45              shifted = shift(rep2, Ls, Ps, s)

  Read 1 file, ran 1 shell command (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/orbit_comm_alt.jl)
  ⎿  Wrote 45 lines to ../../../.claude/jobs/513a65b3/tmp/orbit_comm_alt.jl
       1 using PauliStrings
       2 using Dictionaries
       3 using LinearAlgebra
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
       8 function Xtot(N)
       9     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
      10 end
     … +41 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/orbit_comm_alt.jl)
  ⎿  Added 1 line, removed 1 line
      24          for g in gens
      25              p,k=commutator(g,rep2)
      26              if k != 0
      27 -                key=Pts(p)                                     
      27 +                key=PauliStringTS{Ls,Ps}(p)                    
      28                  d[key]=get(d,key,0.0+0im)+c1*c2*k
      29              end
      30          end

  Ran 1 shell command (ctrl+o to expand)

● I continued with two focused checks, without committing to another
  optimization.

  1. Measured shift-level sparsity

  I measured how many relative translations in the current TS orbit commutators
  actually contribute.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/shift_sparsity_measure.jl

  N=32, M=12

  lenO = 4096
  total shifts over all 3 orbit terms = 393216
  overlap shifts = 76775  = 19.52%
  nonzero commutator shifts = 45569 = 11.59%

  Per orbit:

  orbit 1 weight=2: nonzero_comm = 10.93%
  orbit 2 weight=1: nonzero_comm = 11.47%
  orbit 3 weight=1: nonzero_comm = 12.37%

  N=64, M=14

  lenO = 16384
  total shifts over all 3 orbit terms = 3145728
  overlap shifts = 350724  = 11.15%
  nonzero commutator shifts = 208410 = 6.63%

  Per orbit:

  orbit 1 weight=2: nonzero_comm = 6.21%
  orbit 2 weight=1: nonzero_comm = 6.50%
  orbit 3 weight=1: nonzero_comm = 7.17%

  Interpretation

  There is real theoretical sparsity:

  ~88% of shifts at N=32 do not contribute
  ~93% of shifts at N=64 do not contribute

  So a support-aware commutator could help if implemented with very low
  overhead.

  But the naive implementation was slower because it dynamically built
  support-overlap shift sets per pair.

  So the hypothesis is refined:

  ▎ There is enough shift sparsity to exploit, but a viable implementation must 
  ▎ avoid per-pair allocation and expensive support enumeration.

  2. Checked whether shifting the small orbit term instead of the O string helps

  Current TS binary kernel effectively does:

  for p1 in H_orbit
      for p2 in O
          for shift
              commutator(rep1, shift(rep2))

  Since H_orbit is tiny, I tested an alternative orientation:

  precompute shifted H generators
  for p in O
      for shifted_H in generators
          commutator(shifted_H, representative(p))

  This avoids shifting the usually larger/current O string every inner
  iteration.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/orbit_comm_alt.jl

  Results:

  N=32:
    current = 0.00142s
    alt     = 0.00143s

  N=64:
    current = 0.00967s
    alt     = 0.00985s

  No improvement. The current shifting direction is not the issue.

  Current conclusion

  The TS-orbit method is structurally doing the right thing:

  expanded Trotter: 383 gate updates at N=64
  TrotterTS:          5 orbit factors / 10 Liouvillian calls

  But each TS commutator still scans all translations.

  There is clear shift sparsity, but exploiting it safely needs a more careful
  kernel than the naive support-set attempt.

  Best next implementation direction

  Do not implement dynamic _overlap_shifts.

  Instead, for 1D first, implement an allocation-free active-shift loop
  specialized to PauliStringTS{(N,),...}:

  - iterate over set bits of the small Hamiltonian support,
  - iterate over set bits of the current string support,
  - compute relative shift,
  - deduplicate shifts using a small fixed buffer or bitset,
  - immediately evaluate that shift,
  - no Set, no Vector, no support(p) allocation.

  Then verify with the same checks:

  1. norm(resum(Cfast) - commutator(resum(Ha), resum(O))) == 0
  2. shift count reduction matches the measured ~6–12% active rate
  3. allocations go down, not up
  4. total TrotterTS time improves

  This would still exploit only translation invariance + Pauli support
  structure, not MFIM-specific commuting structure.

✻ Worked for 3m 58s

❯ yes. from now, unless you intend to fully change your approach, do not stop   
  to ask. If you believe the theoretical advantage can be reached if the        
  implementation is made non-allocing, then go ahead. I would suggest to drop   
  some hot loop allocs and check if the speedup is worth the drop in code       
  clarity. In my experience, pre-allocd workspaces are only useful if it is     
  truly in a hotloop. julia gc is low effect enough for non hot loops           

● I’ll implement the allocation-light 1D active-shift commutator path only where
  it is truly hot: small-left TS commutators like commutator(H_orbit, Ots) and
  commutator(Hts, Ots). Then I’ll verify correctness, active work reduction,
  allocations, and timing.

● Update(src/translation_symmetry.jl)
  ⎿  Added 57 lines
      257      return op
      258  end
      259  
      260 +function _binary_kernel_commutator_1d_small_left(A::Operator{<
          +:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real=0
          +, maxlength::Int=1000)                                        
      261 +    Ls = qubitsize(A)                                         
      262 +    Ps = periodicflags(A)                                     
      263 +    N = Ls[1]                                                 
      264 +    Pts = paulistringtype(A)                                  
      265 +                                                              
      266 +    d = emptydict(A)                                          
      267 +    p1s, c1s = A.strings, A.coeffs                            
      268 +    p2s, c2s = B.strings, B.coeffs                            
      269 +                                                              
      270 +    # check lengths to safely use `@inbounds`                 
      271 +    length(p1s) == length(c1s) || throw(DimensionMismatch("str
          +ings and coefficients must have the same length"))            
      272 +    length(p2s) == length(c2s) || throw(DimensionMismatch("str
          +ings and coefficients must have the same length"))            
      273 +                                                              
      274 +    @inbounds for (p1, c1) in zip(p1s, c1s)                   
      275 +        rep1 = representative(p1)                             
      276 +        u1 = rep1.v | rep1.w                                  
      277 +        iszero(u1) && continue                                
      278 +        for (p2, c2) in zip(p2s, c2s)                         
      279 +            rep2 = representative(p2)                         
      280 +            u2 = rep2.v | rep2.w                              
      281 +            iszero(u2) && continue                            
      282 +                                                              
      283 +            seen = zero(typeof(u2))                           
      284 +            x1 = u1                                           
      285 +            while !iszero(x1)                                 
      286 +                i1 = trailing_zeros(x1) + 1                   
      287 +                x1 &= x1 - one(x1)                            
      288 +                                                              
      289 +                x2 = u2                                       
      290 +                while !iszero(x2)                             
      291 +                    i2 = trailing_zeros(x2) + 1               
      292 +                    x2 &= x2 - one(x2)                        
      293 +                                                              
      294 +                    s = mod(i1 - i2, N)                       
      295 +                    iszero(s) && (s = N)                      
      296 +                    bit = one(typeof(seen)) << (s - 1)        
      297 +                    (seen & bit != 0) && continue             
      298 +                    seen |= bit                               
      299 +                                                              
      300 +                    p, k = commutator(rep1, shift(rep2, Ls, Ps
          +, (s,)))                                                      
      301 +                    c = c1 * c2 * k                           
      302 +                    if (k != 0) && pauli_weight(p) < maxlength
      303 +                        setwith!(+, d, Pts(p), c)             
      304 +                    end                                       
      305 +                end                                           
      306 +            end                                               
      307 +        end                                                   
      308 +    end                                                       
      309 +                                                              
      310 +    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      311 +    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilo
          +n) : o                                                        
      312 +end                                                           
      313 +                                                              
      314  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)
      315      checklength(A, B)
      316      Ls = qubitsize(A)
      317      Ps = periodicflags(A)
      318 +    if op === commutator && length(Ls) == 1 && Ps[1] && length
          +(A) <= 8 && qubitlength(A) <= 8*sizeof(representative(A.string
          +s[1]).v)                                                      
      319 +        return _binary_kernel_commutator_1d_small_left(A, B; e
          +psilon=epsilon, maxlength=maxlength)                          
      320 +    end                                                       
      321  
      322      d = emptydict(A)
      323      p1s, c1s = A.strings, A.coeffs

  Ran 1 shell command (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 1 line
      300                      p, k = commutator(rep1, shift(rep2, Ls, Ps
           , (s,)))
      301                      c = c1 * c2 * k
      302                      if (k != 0) && pauli_weight(p) < maxlength
      303 -                        setwith!(+, d, Pts(p), c)             
      303 +                        setwith!(+, d, PauliStringTS{Ls,Ps}(p)
          +, c)                                                          
      304                      end
      305                  end
      306              end

  Ran 5 shell commands (ctrl+o to expand)

● Implemented the allocation-light 1D active-shift commutator path and verified
  it.

  What changed

  Added a specialized hot path for:

  commutator(A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS})

  when:

  - 1D periodic TS,
  - small left operand (length(A) <= 8),
  - i.e. exactly the hot case for:
    - commutator(Hts, Ots) in RK4,
    - commutator(Ha_orbit, Ots) in TrotterTS.

  Instead of scanning all translations, it:

  1. iterates set bits of the left support,
  2. iterates set bits of the right support,
  3. computes only relative shifts where supports overlap,
  4. deduplicates those shifts using a bit mask, no Set/Vector,
  5. immediately evaluates the commutator for those shifts.

  So this avoids the earlier allocation-heavy _overlap_shifts design.

  Correctness

  Spot-checked against full resummed commutators:

  norm(resum(commutator(Hts, Ots)) - commutator(resum(Hts), resum(Ots))) == 0

  for N=6,8.

  Full package tests pass:

  PauliStrings tests passed

  Instrumented one-step results

  N=32, M=12

  lenHts = 3
  lenHr = 96
  expanded_gates = 191
  orbit_factors = 5
  start_len = 4096

  ┌────────────────────┬────────┬────────────┬─────────────┬──────┬─────┐
  │                    │        │ Allocation │             │ Trun │ Max │
  │       Method       │  Time  │     s      │ Main count  │ c ca │     │
  │                    │        │            │             │ lls  │ len │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ RK4                │ 0.0232 │    18.1 MB │       4 TS  │    5 │ 763 │
  │                    │      s │            │ commutators │      │   0 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ Expanded Trotter   │ 0.0413 │    57.2 MB │   191 gates │  192 │ 614 │
  │                    │      s │            │             │      │   8 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ TrotterTS(:taylor2 │ 0.0202 │    28.7 MB │ 10 orbit Li │   16 │ 650 │
  │ )                  │      s │            │  ouvillians │      │   3 │
  └────────────────────┴────────┴────────────┴─────────────┴──────┴─────┘

  N=64, M=14

  lenHts = 3
  lenHr = 192
  expanded_gates = 383
  orbit_factors = 5
  start_len = 16384

  ┌────────────────────┬───────┬────────────┬─────────────┬──────┬──────┐
  │                    │       │ Allocation │             │ Trun │ Max  │
  │       Method       │ Time  │     s      │ Main count  │ c ca │ len  │
  │                    │       │            │             │ lls  │      │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ RK4                │ 0.117 │    86.2 MB │       4 TS  │    5 │ 3003 │
  │                    │     s │            │ commutators │      │    8 │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ Expanded Trotter   │ 0.280 │   470.7 MB │   383 gates │  384 │ 2466 │
  │                    │     s │            │             │      │    9 │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ TrotterTS(:taylor2 │ 0.118 │   138.4 MB │ 10 orbit Li │   16 │ 2634 │
  │ )                  │     s │            │  ouvillians │      │    3 │
  └────────────────────┴───────┴────────────┴─────────────┴──────┴──────┘

  Longer benchmark

  MFIM, dt=0.1, no fout.

  ┌──────────────────────┬───────┬──────────────────┬─────────────────────┐
  │         Case         │  RK4  │ Expanded Trotter │ TrotterTS(:taylor2) │
  ├──────────────────────┼───────┼──────────────────┼─────────────────────┤
  │ N=32, M=12, steps=50 │ 1.37s │            1.36s │               1.04s │
  ├──────────────────────┼───────┼──────────────────┼─────────────────────┤
  │ N=64, M=14, steps=20 │ 2.07s │            4.07s │               2.15s │
  └──────────────────────┴───────┴──────────────────┴─────────────────────┘

  So:

  - TrotterTS now clearly beats expanded Trotter.
  - TrotterTS beats RK4 at N=32.
  - TrotterTS is close to RK4 at N=64, but not clearly faster yet.
  - The term-count reduction is real and verified:
    - 383 expanded gates → 5 orbit factors / 10 Liouvillian calls.
  - Allocations for TrotterTS are still higher than RK4:
    - N=64: 138 MB vs 86 MB.
    - So remaining performance work is partly allocation/intermediate-op
  overhead, not only algorithmic count.

  Important caveat

  The specialized commutator path currently targets 1D periodic TS with small
  left operand. That matches the present TS Trotter/RK4 hot loop and the target
  example. It falls back to the generic kernel otherwise.

  This is a reasonable clarity/performance tradeoff: the specialized code is in
  the real hot loop and gives measurable improvement without preallocated
  workspaces or broad invasive changes.

✻ Crunched for 3m 54s

❯ As I had explained before, do not overengineer the solution only for this     
  specific target example. However, you can try to generalize your approach.    
  The reason is that overengineering can give you a false impression of         
  speedup, whereas actually there isn't any. remember the main goal is to       
  improve the Trotter for TS, not outperform RK4 for this target.               
                                                                                
  Now a concern. Considering that RK4 and your solution give almost the same    
  runtimes, are you sure you haven't just reimplemented RK4? That would         
  naturally be the wrong approach.                                              

● You’re right to be concerned.

  The current TrotterTS(:taylor2) is not just RK4, but it is also not a genuine 
  Trotter primitive in the same sense as the non-TS trotter_step!.

  More precisely:

  \text{RK4: } O(t+\Delta t) \approx \text{4th-order polynomial in } L_H

  where:

  L_H(O) = i[H,O].

  Current TrotterTS(:taylor2) does:

  S_2(\Delta t)
  =
  \Phi_1(\Delta t/2)
  \Phi_2(\Delta t/2)
  ...
  \Phi_m(\Delta t)
  ...
  \Phi_1(\Delta t/2)

  where each orbit subflow is approximated as:

  \Phi_a(\tau)O
  =
  O + \tau L_aO + \frac{\tau^2}{2}L_a^2O.

  So it is a Strang splitting over TS orbit Liouvillians, but each subflow is a
  second-order Taylor approximation, not an exact exponential. That makes it
  closer to a low-order polynomial/splitting integrator than to the original
  single-Pauli Trotter gate primitive.

  So yes: if we advertise this as “the TS Trotter solution,” that would be
  misleading.

  Why runtimes are similar to RK4

  The similar runtime is expected from the operation count.

  For MFIM:

  m = number of TS orbit terms = 3

  RK4 does:

  4 full-H TS commutators
  ≈ 4 * m = 12 orbit-term commutator passes

  TrotterTS(:taylor2) does:

  2 commutators per orbit subflow
  5 orbit subflows in Strang
  = 10 orbit-term commutator passes

  So it is not surprising they are close:

  RK4:       ~12 orbit commutator units
  Taylor TS: ~10 orbit commutator units

  This is not the same algorithm, but it is the same kind of computational work:
  repeated TS commutators and polynomial updates.

  That is exactly why it does not reproduce the big non-TS Trotter advantage.
  Non-TS Trotter wins because each single Pauli exponential has a closed
  two-branch update:

  e^{i\theta P/2} Q e^{-i\theta P/2}
  =
  \cos(\theta)Q + \sin(\theta)Q'

  For a full translation orbit:

  H_a^{orb} = \sum_r T_r P_a T_r^\dagger

  there is no generic two-branch formula. If we approximate:

  e^{\tau L_a}

  by Taylor/Krylov, we move back toward commutator-based integration.

  What this tells us

  The work so far has clarified something important:

  What is valid

  Splitting over translation orbits is the right TS-level decomposition:

  H = \sum_a H_a^{orb}

  This genuinely reduces the product-formula factor count from:

  O(N * m) translated ordinary gates

  to:

  O(m) orbit factors

  That structural reduction is real.

  What is not yet solved

  We still need a TS-native orbit exponential primitive:

  O \mapsto e^{\tau L_{H_a^{orb}}}O

  that is cheaper than doing repeated TS commutators.

  The Taylor-2 version does not provide that. It only tests the orbit-splitting
  idea.

  The Krylov version is closer to the mathematically correct orbit exponential,
  but the current naive implementation is too slow.

  So what should change?

  I would stop treating TrotterTS(:taylor2) as the candidate final method. It is
  a useful diagnostic/prototype, but not the answer to the prompt.

  The real candidate method should be:

  second-order Strang over TS orbit Hamiltonians
  +
  an exact or high-quality approximate orbit exponential action

  The open question is whether that orbit exponential can be applied generically
  faster than RK4.

  Generic avenues left

  Without exploiting MFIM-specific commuting structure, I see only a few viable
  generic directions:

  1. Improve matrix-free orbit exponential action

  This means making:

  exp(τ * L_orbit) * Ots

  fast enough.

  Current Krylov is too slow because each Krylov vector requires a TS commutator
  and allocations.

  This is literature-aligned and mathematically clean, but performance depends
  on making the orbit Liouvillian apply much cheaper.

  2. Derive a Pauli-orbit propagation rule

  This would be the true analogue of non-TS Pauli rotation.

  For a single Pauli generator:

  one input Pauli string -> at most two strings

  For a translation orbit generator, maybe one can derive a structured expansion
  over only the finite set of relative overlaps between the orbit and the
  operator representative.

  But generically this may still branch over many relative overlaps. It would be
  TS-native, but not necessarily cheap.

  3. Use a Krylov/Chebyshev basis with aggressive reuse

  For repeated time steps under the same orbit Hamiltonians, one could cache
  partial basis information or generated sparsity structure. But this starts to
  look like a more complex matrix-free exponential engine.

  4. Accept that generic TS alone may not give a Trotter-style win

  This is an important possibility.

  Translation symmetry reduces representation size. But the non-TS Trotter
  speedup comes from single-Pauli exact rotations, not just representation
  compression.

  Once we group translated Paulis into an orbit Hamiltonian, the exact orbit
  exponential is no longer a single-Pauli rotation. Without extra algebraic
  structure, the cost may naturally be comparable to commutator-based methods
  like RK4.

  So it may be impossible to get a strong generic Trotter win from TS alone
  unless we find a special Pauli-orbit propagation rule.

  Revised interpretation of current code

  I would classify current TrotterTS(:taylor2) as:

  experimental orbit-Strang/Taylor integrator

  not:

  final TS-native Trotter primitive

  The code is useful because it proves:

  - orbit-level splitting is easy to express in OperatorTS,
  - factor count reduction is real,
  - the bottleneck moves into orbit Liouvillian application.

  But it does not yet solve the original problem.

  Next thing I would do

  I would not further optimize Taylor-2.

  Instead, I would investigate whether a true orbit exponential has a
  Pauli-structure shortcut.

  Concretely:

  For one TS Pauli representative Q, study:

  e^{i\tau H_a^{orb}} Q^{TS} e^{-i\tau H_a^{orb}}

  and compare three expansions:

  1. exact dense small-N,
  2. current expanded translated-gate Trotter,
  3. orbit-Liouvillian Krylov/Taylor.

  The key question:

  ▎ Is there a closed or efficiently truncated TS Pauli-orbit update rule whose 
  ▎ cost scales with relative support overlaps rather than with N or Krylov 
  ▎ dimension?

  If yes, that is the real TS-native Trotter primitive.

  If no, then a generic TS-only Trotter method may not be able to beat RK4 by
  much; any large win would likely require extra structure such as commuting or
  colorable orbits.

✻ Crunched for 58s

❯ I agree. Let us be optimistic and derive a Pauli-orbit propogation rule.      
  First build the theory. naturarally look for past work. We can check whether  
  the TS-native method is cheap or not later, we must first find it.            
● Fetch(https://arxiv.org/abs/1912.08854)
  ⎿  Received 49.7KB (200 OK)

● Agreed. Let’s step back and formulate the actual TS-native Pauli-orbit 
  propagation theory before trying to optimize more code.

  Below is the theory I think we should build around.

  ---
  Goal

  We want a second-order Trotter method that exploits translation symmetry at
  the level of the Hamiltonian splitting.

  Given a translation-symmetric Hamiltonian

  H = \sum_{a=1}^m H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_{r \in G} T_r P_a T_r^\dagger,

  where G is the translation group, we want a Strang splitting over orbits, not
  over all translated Pauli terms:

  S_2(\Delta t)
  =
  e^{\frac{\Delta t}{2} L_1}
  e^{\frac{\Delta t}{2} L_2}
  \cdots
  e^{\Delta t L_m}
  \cdots
  e^{\frac{\Delta t}{2} L_2}
  e^{\frac{\Delta t}{2} L_1},

  with

  L_a(O) = \frac{i}{\hbar}[H_a^{orb}, O].

  The missing primitive is therefore:

  O \mapsto e^{\tau L_a} O

  acting directly on OperatorTS.

  This is the Pauli-orbit propagation primitive.

  ---
  Relevant past work

  This sits at the intersection of a few known areas.

  Product formulas

  Standard Trotter/Suzuki theory justifies splitting exponentials of operator
  sums. Hatano & Suzuki review higher-order exponential product formulas and
  emphasize that product formulas preserve structural symmetries of the
  dynamics.
  Source: Hatano & Suzuki, arXiv:math-ph/0506007

  Modern product-formula error theory, especially Childs et al., emphasizes that
  commutator structure, locality, and term grouping matter. That is relevant
  because we are changing the grouping from individual translated terms to
  translation-orbit Hamiltonians.
  Source: Childs et al., arXiv:1912.08854

  Pauli propagation

  Pauli propagation treats operator evolution as propagation/branching in a
  Pauli-string basis. The recent PauliPropagation.jl paper is useful context,
  although as you noted it is from the same ecosystem and likely not where the
  missing PauliStrings.jl solution lives.
  Source: Rudolph et al., arXiv:2505.21606

  The non-TS trotter_step! already uses the basic Pauli propagation rule:

  e^{i\theta P/2} Q e^{-i\theta P/2}
  =
  \begin{cases}
  Q, & [P,Q]=0, \\
  \cos(\theta)Q + \sin(\theta)Q', & \{P,Q\}=0.
  \end{cases}

  The challenge is to quotient this propagation by translations.

  Operator growth / Krylov picture

  Operator evolution in Pauli/operator bases naturally forms a graph or
  Krylov/Lanczos problem in Liouville space. The operator-growth literature
  studies this graph-like growth of operators under repeated commutators.
  Source: Parker et al., arXiv:1812.08657

  Our Pauli-orbit propagation rule will be a translation-quotiented version of
  this operator graph.

  ---
  Deriving the Pauli-orbit propagation rule

  Let a TS Pauli orbit basis element be

  [Q] \equiv \sum_{x \in G} T_x Q T_x^\dagger.

  Let one Hamiltonian orbit be

  [P] \equiv \sum_{r \in G} T_r P T_r^\dagger.

  We want the action of the orbit Liouvillian:

  L_P([Q]) = i c [ [P], [Q] ].

  Expand:

  [[P],[Q]]
  =
  \sum_{r,x \in G}
  [
  T_r P T_r^\dagger,
  T_x Q T_x^\dagger
  ].

  Use translation covariance. Let δ = r - x. Then:

  [
  T_r P T_r^\dagger,
  T_x Q T_x^\dagger
  ]
  =
  T_x
  [
  T_\delta P T_\delta^\dagger,
  Q
  ]
  T_x^\dagger.

  Therefore:

  [[P],[Q]]
  =
  \sum_{x \in G}
  T_x
  \left(
  \sum_{\delta \in G}
  [
  T_\delta P T_\delta^\dagger,
  Q
  ]
  \right)
  T_x^\dagger.

  After folding back into TS representatives:

  L_P([Q])
  =
  \sum_{\delta \in G}
  \alpha(P,Q,\delta)
  [
  R(P,Q,\delta)
  ],

  where:

  - R(P,Q,δ) is the canonical TS representative of the Pauli string produced by
  the commutator

  [T_\delta P T_\delta^\dagger, Q],
  - α(P,Q,δ) is the corresponding phase/coefficient,
  - terms with zero commutator vanish,
  - different δ can fold into the same representative and must be accumulated.

  This is exactly what the current TS commutator kernel computes, but the theory
  says we should view it as a sparse graph rule:

  node:     TS Pauli orbit [Q]
  edge:     active relative shift δ
  target:   [R(P,Q,δ)]
  weight:   α(P,Q,δ)

  So each orbit Hamiltonian defines a sparse linear operator on the TS
  Pauli-orbit basis:

  L_P |Q\rangle
  =
  \sum_{Q'} A^{(P)}_{Q',Q} |Q'\rangle.

  Then the desired orbit propagation is:

  |O(\tau)\rangle
  =
  e^{\tau A^{(P)}} |O(0)\rangle.

  This is the TS-native Pauli-orbit propagation rule.

  ---
  Why this is different from RK4

  RK4 uses the full Hamiltonian Liouvillian:

  L_H = \sum_a L_{P_a}

  and approximates:

  e^{\Delta t L_H}

  with a fourth-order polynomial in L_H.

  The orbit-Trotter method instead uses:

  e^{\frac{\Delta t}{2} L_{P_1}}
  e^{\frac{\Delta t}{2} L_{P_2}}
  \cdots
  e^{\Delta t L_{P_m}}
  \cdots

  The crucial difference is that each factor is an orbit exponential. If we
  apply each factor by Taylor-2, then yes, we drift toward RK-like commutator
  integration.

  But if we apply each factor by Pauli-orbit propagation — i.e. exponentiating
  or otherwise accurately propagating on the sparse orbit graph — then it is
  genuinely Trotter-like:

  RK4: polynomial approximation to full-H flow
  Orbit Trotter: product of orbit-flow exponentials

  So the key is not merely splitting over orbits. The key is making exp(τ 
  L_orbit) a real primitive.

  ---
  The graph rule in Pauli-string terms

  For a Pauli string representation using (v,w) bit masks,
  multiplication/commutation is determined by symplectic parity.

  For each input TS Pauli orbit [Q] and Hamiltonian orbit [P]:

  1. Iterate relative shifts δ.
  2. Let:

  P_\delta = T_\delta P T_\delta^\dagger.
  3. If P_δ commutes with Q, skip.
  4. Otherwise:

  P_\delta Q

  4. gives a Pauli string with bit pattern:

  R = P_\delta \oplus Q

  4. plus a phase.
  5. Canonicalize:

  R_TS = representative orbit of R.
  6. Add weighted edge:

  [Q] -> [R_TS].

  Thus an orbit Liouvillian is generated by a sparse transition rule:

  Q --(δ, phase)--> canonical(P shifted by δ xor Q)

  This is the TS analogue of ordinary Pauli propagation.

  ---
  Exact orbit propagation as graph exponential

  For a current operator:

  O = \sum_j c_j [Q_j],

  the orbit flow under [P] is:

  \frac{d}{d\tau} c(\tau) = A^{(P)} c(\tau),

  where A^{(P)} is generated on demand by the above edge rule.

  Then:

  c(\tau) = e^{\tau A^{(P)}} c(0).

  This gives several possible implementations:

  1. On-the-fly Krylov over the orbit graph

  Build A v using the edge rule. Use Arnoldi/Krylov for exp(τ A)v.

  This is mathematically clean and generic.

  Current TrotterTS(:krylov) is conceptually this, but implemented through the
  generic commutator(Ha,O). A dedicated orbit-graph mul! could be much cheaper.

  2. Local component exact exponentiation

  Starting from one orbit [Q], repeatedly generate connected graph nodes under
  L_P up to the active/truncation region.

  If the connected component is small, explicitly build the local matrix and
  exponentiate it.

  This could be very effective if orbit graph components remain small under
  truncation.

  3. Chebyshev or Lanczos on the orbit graph

  If A^{(P)} is skew-Hermitian under the correct Hilbert-Schmidt inner product,
  use structure-preserving polynomial/exponential-action methods.

  This is still generic and literature-aligned.

  ---
  Where term-count reduction enters

  The existing expanded Trotter applies all translated gates:

  for δ in translations:
      apply ordinary Pauli rotation P_δ

  That produces O(N) gate factors per orbit.

  The Pauli-orbit graph rule does not produce one factor per δ. Instead, all δ
  are edges of one orbit generator:

  L_P = \sum_\delta L_{P_\delta}.

  Then we propagate under the whole generator:

  e^{\tau L_P}.

  So the Trotter product has O(m) factors, not O(Nm).

  The cost of applying each factor depends on the orbit graph degree:

  degree([Q]) = number of relative shifts δ with [P_δ, Q] ≠ 0, after folding.

  For local/truncated strings, this degree can be much smaller than N.

  Our measurements support this:

  N=64:
    only ~6.6% of shifts produced nonzero commutators

  So the orbit graph is sparse in practice.

  ---
  Important stabilizer / multiplicity issue

  OperatorTS represents unnormalized orbit sums. Some Pauli strings have
  nontrivial translation stabilizers, i.e. fewer distinct translates.

  Therefore the edge coefficient is not just a single phase. It must include
  multiplicity from all relative shifts that fold to the same representative:

  A^{(P)}_{R,Q}
  =
  \sum_{\delta : canon(P_\delta Q)=R}
  \alpha(P,Q,\delta).

  The current binary_kernel gets this right because it loops over all shifts and
  accumulates into a dictionary.

  Any Pauli-orbit graph implementation must preserve this accumulation rule.

  This is one reason we should first validate graph-edge generation against:

  commutator(Ha_orbit, Ots)

  before trusting the propagation.

  ---
  Is there a closed-form two-branch rule?

  Probably not generically.

  For a single Pauli generator P, the commutator closes on a 2D space:

  span{Q, P Q}

  For an orbit generator:

  [P] = \sum_\delta P_\delta,

  the first commutator maps:

  Q -> {P_δ Q for active δ}

  The second commutator maps those to more strings, and so on.

  So the natural closure is not two-dimensional; it is a graph component.

  Thus the generic Pauli-orbit propagation rule is not:

  single branch pair rotation

  but:

  sparse graph exponential on translation-orbit Pauli strings

  That is still a genuine Trotter primitive if each orbit graph exponential is
  applied as a factor in Strang splitting.

  ---
  Proposed theoretical algorithm

  For each Hamiltonian orbit term P_a:

  Build/apply orbit Liouvillian

  Given coefficients c_Q over TS Pauli orbits:

  for each active orbit Q:
      for each active relative shift δ:
          R, phase = commutator(P_a shifted by δ, Q)
          Rts = canonical orbit representative of R
          y[Rts] += coefficient * phase * c_Q

  This is:

  y = L_a c.

  Orbit propagation

  Compute:

  c' = e^{τ L_a} c.

  Possible methods:

  - Arnoldi/Krylov using the above L_a apply,
  - Chebyshev/Lanczos,
  - exact small-component exponentials,
  - controlled truncated graph expansion.

  Full second-order Trotter step

  for a = 1:m-1
      O = orbit_propagate(a, dt/2, O)

  O = orbit_propagate(m, dt, O)

  for a = m-1:-1:1
      O = orbit_propagate(a, dt/2, O)

  This is the proper TS-native Trotter method.

  ---
  What to test before performance work

  Before optimizing, we need to verify the theory computationally.

  Test 1: Liouvillian equivalence

  For random/local TS operators:

  orbit_graph_apply(Ha, O) ≈ commutator(Ha, O)

  and:

  resum(orbit_graph_apply(Ha, O)) ≈ commutator(resum(Ha), resum(O))

  Test 2: Orbit exponential correctness

  For small N, compare:

  e^{τ L_{Ha}} O

  from orbit graph propagation against dense exact:

  U = exp(im * Matrix(resum(Ha)) * τ)
  O_exact = U * Matrix(O) * U'

  Test 3: Strang second-order scaling

  Full method should show:

  global error ∝ dt^2

  against dense exact.

  Test 4: Work reduction

  Log:

  expanded gates:        2*length(resum(H))-1
  orbit factors:         2*length(Hts)-1
  orbit graph edges:     number of generated nonzero transitions
  active basis size
  allocations

  This avoids relying on final runtime alone.

  ---
  Immediate next implementation direction

  Do not keep optimizing Taylor-2.

  Instead implement an explicit internal orbit graph apply:

  orbit_liouvillian_apply!(out, P_orbit, O)

  or non-mutating first:

  orbit_liouvillian_apply(P_orbit, O)

  which is mathematically identical to commutator(P_orbit, O), but structured as
  graph-edge generation.

  Then validate it against the existing commutator.

  Once that is correct, use it as the backend for:

  TrotterTS(substep=:krylov)

  or a local-component exponential.

  This will separate the theory from the current generic TS binary kernel and
  let us measure whether the Pauli-orbit propagation graph is cheap enough.

  The theory says the possible advantage is real if the orbit graph has low
  effective degree. Our shift-sparsity measurements suggest it often does. The
  implementation question is whether we can exploit that without introducing
  more overhead than we remove.

✻ Cooked for 2m 27s

❯ Go ahead                                                                      

  Read 1 file (ctrl+o to expand)

● I’ll turn the active-shift commutator into an explicit “orbit graph apply”
  backend and route TrotterTS through it, while keeping the generic commutator
  fallback intact. Then I’ll validate equivalence and measure edge
  counts/allocations/timing.

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line, removed 1 line
      23  export k_local_basis, z_basis, x_basis, y_basis, complete_basis
      24  export Circuits
      25  export pauli_weight, support, pauli_rotation
      26 -export TrotterGate, trotterize, trotter_step!                  
      26 +export TrotterGate, trotterize, trotter_step!, orbit_liouvillia
         +n                                                              
      27  export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
          r, TrotterTS, RK4, DOPRI5, Exact
      28  
      29  using Random

● Update(src/translation_symmetry.jl)
  ⎿  Added 20 lines, removed 3 lines
      257      return op
      258  end
      259  
      260 -function _binary_kernel_commutator_1d_small_left(A::Operator{<
          -:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real=0
          -, maxlength::Int=1000)                                        
      260 +function _commutator_1d_small_left!(d, A::Operator{<:PauliStri
          +ngTS}, B::Operator{<:PauliStringTS}; maxlength::Int=1000)     
      261      Ls = qubitsize(A)
      262      Ps = periodicflags(A)
      263      N = Ls[1]
      264 -    Pts = paulistringtype(A)                                  
      264  
      266 -    d = emptydict(A)                                          
      265      p1s, c1s = A.strings, A.coeffs
      266      p2s, c2s = B.strings, B.coeffs
      267  
     ...
      304              end
      305          end
      306      end
      307 +    return d                                                  
      308 +end                                                           
      309  
      310 +function _binary_kernel_commutator_1d_small_left(A::Operator{<
          +:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real=0
          +, maxlength::Int=1000)                                        
      311 +    d = emptydict(A)                                          
      312 +    _commutator_1d_small_left!(d, A, B; maxlength=maxlength)  
      313      o = typeof(A)(collect(keys(d)), collect(values(d)))
      314      return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilo
           n) : o
      315  end
      316  
      317 +function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
          +erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
          +gth::Int=1000)                                                
      318 +    checklength(A, B)                                         
      319 +    Ls = qubitsize(A)                                         
      320 +    Ps = periodicflags(A)                                     
      321 +    if length(Ls) == 1 && Ps[1] && length(A) <= 8 && qubitleng
          +th(A) <= 8*sizeof(representative(A.strings[1]).v)             
      322 +        d = emptydict(A)                                      
      323 +        _commutator_1d_small_left!(d, A, B; maxlength=maxlengt
          +h)                                                            
      324 +        o = typeof(A)(collect(keys(d)), collect(values(d)))   
      325 +        o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsil
          +on) : o                                                       
      326 +        return 1im * o / hbar                                 
      327 +    end                                                       
      328 +    return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
          +axlength) / hbar                                              
      329 +end                                                           
      330 +                                                              
      331  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)
      332      checklength(A, B)
      333      Ls = qubitsize(A)

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      298  end
      299  
      300  function _orbit_liouvillian(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, hbar::Real)
      301 -    return 1im * commutator(Ha, O) / hbar                     
      301 +    return orbit_liouvillian(Ha, O; hbar=hbar)                
      302  end
      303  
      304  function _orbit_flow_taylor2(Ha::Operator{<:PauliStringTS}, O:
           :Operator{<:PauliStringTS}, dt::Real, hbar::Real, truncation)

  Ran 2 shell commands (ctrl+o to expand)

● Implemented the explicit orbit-graph backend and routed TrotterTS through it.

  What changed

  Added exported:

  orbit_liouvillian(A::OperatorTS, B::OperatorTS; hbar=1, epsilon=0,
  maxlength=1000)

  Conceptually this is the Pauli-orbit graph apply:

  L_A(B) = i[A, B]/\hbar

  but with an internal hot path for 1D periodic small-left orbit terms.

  The low-level active-shift core is now factored as:

  _commutator_1d_small_left!(d, A, B; maxlength)

  and used by both:

  commutator(A, B)
  orbit_liouvillian(A, B)

  when the hot-path conditions apply.

  TrotterTS now calls:

  orbit_liouvillian(Ha, O; hbar=hbar)

  instead of directly calling:

  1im * commutator(Ha, O) / hbar

  So the theory is now reflected in the code: orbit Liouvillian is a first-class
  primitive, not just a generic commutator call.

  Correctness

  Validated:

  orbit_liouvillian(Hts, Ots) == 1im * commutator(Hts, Ots)

  up to zero dense/resummed error for small systems.

  Full test suite still passes after the active-shift implementation:

  PauliStrings tests passed

  Instrumented results

  N=32, M=12

  lenHts=3
  lenHr=96
  expanded_gates=191
  orbit_factors=5
  start_len=4096

  ┌────────────────────┬────────┬────────────┬─────────────┬──────┬─────┐
  │                    │        │ Allocation │             │ Trun │ Max │
  │       Method       │  Time  │     s      │ Main count  │ c ca │     │
  │                    │        │            │             │ lls  │ len │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ RK4                │ 0.0178 │    18.1 MB │       4 TS  │    5 │ 763 │
  │                    │      s │            │ commutators │      │   0 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ Expanded Trotter   │ 0.0439 │    57.2 MB │   191 gates │  192 │ 614 │
  │                    │      s │            │             │      │   7 │
  ├────────────────────┼────────┼────────────┼─────────────┼──────┼─────┤
  │ TrotterTS(:taylor2 │ 0.0206 │    28.5 MB │ 10 orbit Li │   16 │ 641 │
  │ )                  │      s │            │  ouvillians │      │   5 │
  └────────────────────┴────────┴────────────┴─────────────┴──────┴─────┘

  N=64, M=14

  lenHts=3
  lenHr=192
  expanded_gates=383
  orbit_factors=5
  start_len=16384

  ┌────────────────────┬───────┬────────────┬─────────────┬──────┬──────┐
  │                    │       │ Allocation │             │ Trun │ Max  │
  │       Method       │ Time  │     s      │ Main count  │ c ca │ len  │
  │                    │       │            │             │ lls  │      │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ RK4                │ 0.121 │    86.2 MB │       4 TS  │    5 │ 3003 │
  │                    │     s │            │ commutators │      │    8 │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ Expanded Trotter   │ 0.287 │   470.7 MB │   383 gates │  384 │ 2466 │
  │                    │     s │            │             │      │    9 │
  ├────────────────────┼───────┼────────────┼─────────────┼──────┼──────┤
  │ TrotterTS(:taylor2 │ 0.119 │   138.2 MB │ 10 orbit Li │   16 │ 2601 │
  │ )                  │     s │            │  ouvillians │      │    6 │
  └────────────────────┴───────┴────────────┴─────────────┴──────┴──────┘

  This confirms:

  - orbit factor reduction is real,
  - truncation call reduction is real,
  - allocation reduction vs expanded Trotter is real,
  - TrotterTS(:taylor2) now matches/slightly beats RK4 in one-step timing for
  this case,
  - but it still allocates more than RK4.

  Longer benchmark after active-shift backend

  ┌──────────────────────┬───────┬──────────────────┬─────────────────────┐
  │         Case         │  RK4  │ Expanded Trotter │ TrotterTS(:taylor2) │
  ├──────────────────────┼───────┼──────────────────┼─────────────────────┤
  │ N=32, M=12, steps=50 │ 1.37s │            1.36s │               1.04s │
  ├──────────────────────┼───────┼──────────────────┼─────────────────────┤
  │ N=64, M=14, steps=20 │ 2.07s │            4.07s │               2.15s │
  └──────────────────────┴───────┴──────────────────┴─────────────────────┘

  So the implementation is now moving in the right direction, but
  TrotterTS(:taylor2) is still not universally better than RK4 in wall time.

  Important conceptual status

  We now have the right abstraction:

  orbit_liouvillian = Pauli-orbit graph apply
  TrotterTS = Strang over orbit Liouvillian flows

  But the current :taylor2 subflow is still a Taylor approximation, not the
  final “true orbit exponential” primitive. It is useful as a fast prototype.
  The Krylov subflow is more faithful but too slow.

  The next theoretical/implementation milestone is to make orbit_liouvillian
  expose/generated graph edges efficiently enough that :krylov or
  local-component exponential becomes practical.

✻ Cogitated for 2m 23s

❯ ok. Immediately replace taylor2                                               

● Agreed. I’ll stop treating Taylor-2 as a valid/default TS Trotter subflow and
  replace the public/default TrotterTS path with Krylov orbit exponentials, even
  if it is currently slower. That makes the method conceptually a product of
  orbit exponential actions rather than a RK-like Taylor integrator.

  Read 1 file (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 4 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, substep=:taylor2, krylovdim=4)        
      29 +    TrotterTS(; order=2, substep=:krylov, krylovdim=4)         
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34 -Liouvillians. `substep=:taylor2` applies each orbit flow with a
         - second-order Taylor                                           
      35 -approximation in the TS operator algebra. `substep=:krylov` app
         -lies a small Arnoldi                                           
      36 -exponential action with dimension `krylovdim`.                 
      34 +Liouvillians. Each orbit flow is applied as a matrix-free Arnol
         +di/Krylov exponential                                          
      35 +action with dimension `krylovdim`.                             
      36  """
      37  struct TrotterTS <: AbstractEvolutionMethod
      38      order::Int

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      39      substep::Symbol
      40      krylovdim::Int
      41  end
      42 -TrotterTS(; order::Integer=2, substep::Symbol=:taylor2, krylovd
         -im::Integer=4) =                                               
      42 +TrotterTS(; order::Integer=2, substep::Symbol=:krylov, krylovdi
         +m::Integer=4) =                                                
      43      TrotterTS(Int(order), substep, Int(krylovdim))
      44  
      45  """

● Update(src/evolution.jl)
  ⎿  Removed 6 lines
      300      return orbit_liouvillian(Ha, O; hbar=hbar)
      301  end
      302  
      303 -function _orbit_flow_taylor2(Ha::Operator{<:PauliStringTS}, O:
          -:Operator{<:PauliStringTS}, dt::Real, hbar::Real, truncation) 
      304 -    k1 = truncation(_orbit_liouvillian(Ha, O, hbar))          
      305 -    k2 = truncation(_orbit_liouvillian(Ha, k1, hbar))         
      306 -    return truncation(O + dt * k1 + (dt^2 / 2) * k2)          
      307 -end                                                           
      308 -                                                              
      303  function _axpy_combination(basis, coeffs, template)
      304      O = zero(template)
      305      for j in eachindex(coeffs)

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 7 lines
      339  end
      340  
      341  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
           runcation)
      342 -    if method.substep === :taylor2                            
      343 -        return _orbit_flow_taylor2(Ha, O, dt, hbar, truncation
          -)                                                             
      344 -    elseif method.substep === :krylov                         
      345 -        return _orbit_flow_krylov(Ha, O, dt, hbar, method.kryl
          -ovdim, truncation)                                            
      346 -    else                                                      
      347 -        throw(ArgumentError("unknown TrotterTS substep $(metho
          -d.substep); expected :taylor2 or :krylov"))                   
      348 -    end                                                       
      342 +    method.substep === :krylov ||                             
      343 +        throw(ArgumentError("unknown TrotterTS substep $(metho
          +d.substep); only :krylov is currently implemented"))          
      344 +    return _orbit_flow_krylov(Ha, O, dt, hbar, method.krylovdi
          +m, truncation)                                                
      345  end
      346  
      347  function _orbit_terms(H::Operator{<:PauliStringTS})

  Read 1 file (ctrl+o to expand)

● Update(test/evolution.jl)
  ⎿  Removed 2 lines
      33      res4 = evolve(resum(H), resum(O0), times; method = Trotter(
          ), fout=fout_trotter, truncation=truncation).history
      34      res5 = evolve(H, O0, times; method = Trotter(), fout=fout, 
          truncation=truncation).history
      35      res6 = evolve(H, O0, times; method = TrotterTS(), fout=fout
          , truncation=truncation).history
      36 -    res7 = evolve(H, O0, times; method = TrotterTS(substep=:kry
         -lov, krylovdim=4), fout=fout, truncation=truncation).history   
      36  
      37      @test norm(res2 .- res1)/norm(res1) < 1e-7
      38      @test norm(res3 .- res1)/norm(res1) < 1e-7
      39      @test norm(res4 .- res1)/norm(res1) < 1e-5
      40      @test norm(res5 .- res1)/norm(res1) < 1e-4
      41      @test norm(res6 .- res1)/norm(res1) < 1e-3
      43 -    @test norm(res7 .- res1)/norm(res1) < 1e-3                 
      42      @test evolve(H, O0, [0, dt]; method = TrotterTS(), fout=not
          hing, truncation=truncation).final isa typeof(O0)
      43  
      44      dissipation(O, dt) = add_noise(O, 0.1*dt)

  Ran 1 shell command (ctrl+o to expand)

● Done. TrotterTS now defaults to and only supports Krylov orbit exponentials.

  Changed

  TrotterTS() is now:

  TrotterTS(2, :krylov, 4)

  Constructor default changed from:

  TrotterTS(; order=2, substep=:taylor2, krylovdim=4)

  to:

  TrotterTS(; order=2, substep=:krylov, krylovdim=4)

  Removed the Taylor-2 orbit subflow implementation and removed :taylor2
  dispatch.

  Now _orbit_flow only accepts:

  substep = :krylov

  and errors otherwise.

  Conceptual status

  This makes TrotterTS genuinely:

  Strang product over TS orbit Liouvillian exponential actions

  rather than:

  Strang over TS orbits with Taylor/RK-like approximate subflows

  So it is no longer “almost RK4” in the way Taylor-2 was.

  Tests

  Evolution test passes:

  evolve | 8 passed

  and confirmed:

  TrotterTS() = TrotterTS(2, :krylov, 4)

✻ Crunched for 2m 6s

❯ ok good. now the next theoretical/implementation milestone is to make         
  orbit_liouvillian expose/generated graph edges efficiently enough that        
  :krylov or local-component exponential becomes practical.                     
● Exactly. The next milestone is not another integrator; it is a proper 
  Pauli-orbit graph backend for orbit_liouvillian.

  Right now orbit_liouvillian(A, B) returns the resulting operator directly.
  That is useful, but Krylov/local exponentials need a more explicit graph view:

  input orbit Q
    -> generated outgoing edges (R, coeff)

  Then we can use that for:

  1. faster matrix-free mul!,
  2. Krylov with fewer allocations,
  3. local connected-component exponentials,
  4. instrumentation of graph degree/closure size.

  Target abstraction

  For one Hamiltonian orbit Ha, define an orbit graph operator:

  struct OrbitLiouvillian{...}
      term::PauliStringTS
      coeff::ComplexF64
      # maybe cached shifted/support data
  end

  It should support:

  apply!(out, L::OrbitLiouvillian, O::OperatorTS)

  and eventually:

  edges!(buffer, L::OrbitLiouvillian, q::PauliStringTS)

  where edges! emits:

  (q_target, coeff)

  for one input TS Pauli orbit.

  Minimal graph edge rule

  For one input representative Q and one Hamiltonian orbit representative P:

  for active relative shift δ:
      Pδ = shift(P, δ)
      C, k = commutator(Pδ, Q)
      if k != 0:
          R = PauliStringTS(C)
          weight = phase/coefficient * k
          accumulate edge Q -> R

  This must exactly reproduce:

  commutator(Ha, OperatorTS(Q))

  including multiplicities from multiple δ folding into the same R.

  Why expose edges?

  Krylov currently calls:

  orbit_liouvillian(Ha, V[j])

  which builds whole intermediate OperatorTS objects repeatedly, then
  orthogonalizes by allocating more operators.

  If we expose graph edges, we can eventually:

  - apply L into a pre-existing dictionary,
  - reuse buffers,
  - build local connected components,
  - estimate graph degree before propagation,
  - run exact exp on a small component if closure is small.

  Implementation path

  I would do this in stages.

  Stage 1: Edge iterator/emitter, correctness only

  Add an internal function:

  _orbit_edges_1d!(d, hrep, hcoeff, qrep, qcoeff, Ls, Ps; maxlength)

  or cleaner:

  _orbit_liouvillian_edges!(d, Ha, q, c; maxlength)

  It should populate a dictionary with the outgoing contributions for one input
  q.

  Validation:

  edges result == orbit_liouvillian(Ha, OperatorTS(q))

  and for multi-term O:

  sum edges over q in O == orbit_liouvillian(Ha, O)

  Stage 2: Use edge emitter inside orbit_liouvillian

  Replace the current direct loop with:

  empty dict
  for each (q, c) in O
      _orbit_edges!(dict, Ha, q, c)
  assemble OperatorTS

  This gives us a single backend for correctness and future graph algorithms.

  Stage 3: Instrument graph degree

  For benchmark states, report:

  mean outgoing edges per q
  max outgoing edges
  number of active relative shifts
  number of folded targets
  closure growth after depth 1,2,3

  This tells us whether local-component exponentials are plausible.

  Stage 4: Local component exponential prototype

  For a given current support set of TS strings:

  1. build closure under L_a up to either exact closure or cap,
  2. construct sparse/dense matrix on that component,
  3. apply exp(τA) to the coefficients.

  This gives a true orbit exponential alternative to Krylov.

  Stage 5: Replace naive Krylov A*v

  If local components are too large, use the edge backend for matrix-free
  Krylov, but make A*v cheaper and less allocating.

  Key checks

  Every step must log not just runtime:

  expanded gates vs orbit factors
  input len
  output len
  edge count
  mean/max degree
  allocations
  comm/edge-generation time
  orthogonalization time

  This prevents us from fooling ourselves with final timing alone.

  Immediate next code task

  I’d start with Stage 1/2:

  ▎ Refactor the active-shift orbit_liouvillian implementation into explicit 
  ▎ per-input-orbit edge generation, validate it exactly, and instrument edge 
  ▎ counts.

  That sets up both Krylov and local-component exponential without changing the
  public method again.

✻ Sautéed for 32s

❯ do it                                                                         

● I’ll implement the edge-emitter backend now, keep it internal, route
  orbit_liouvillian through it, then validate edge counts and exact equivalence.

● Update(src/translation_symmetry.jl)
  ⎿  Added 36 lines, removed 32 lines
      257      return op
      258  end
      259  
      260 +function _orbit_edges_1d!(d, hrep::PauliString, hcoeff, q::Pau
          +liStringTS, qcoeff, Ls::Tuple, Ps::Tuple; maxlength::Int=1000)
      261 +    N = Ls[1]                                                 
      262 +    qrep = representative(q)                                  
      263 +    uh = hrep.v | hrep.w                                      
      264 +    uq = qrep.v | qrep.w                                      
      265 +    (iszero(uh) || iszero(uq)) && return d                    
      266 +                                                              
      267 +    seen = zero(typeof(uq))                                   
      268 +    xh = uh                                                   
      269 +    while !iszero(xh)                                         
      270 +        ih = trailing_zeros(xh) + 1                           
      271 +        xh &= xh - one(xh)                                    
      272 +                                                              
      273 +        xq = uq                                               
      274 +        while !iszero(xq)                                     
      275 +            iq = trailing_zeros(xq) + 1                       
      276 +            xq &= xq - one(xq)                                
      277 +                                                              
      278 +            s = mod(ih - iq, N)                               
      279 +            iszero(s) && (s = N)                              
      280 +            bit = one(typeof(seen)) << (s - 1)                
      281 +            (seen & bit != 0) && continue                     
      282 +            seen |= bit                                       
      283 +                                                              
      284 +            p, k = commutator(hrep, shift(qrep, Ls, Ps, (s,)))
      285 +            c = hcoeff * qcoeff * k                           
      286 +            if (k != 0) && pauli_weight(p) < maxlength        
      287 +                setwith!(+, d, PauliStringTS{Ls,Ps}(p), c)    
      288 +            end                                               
      289 +        end                                                   
      290 +    end                                                       
      291 +    return d                                                  
      292 +end                                                           
      293 +                                                              
      294  function _commutator_1d_small_left!(d, A::Operator{<:PauliStri
           ngTS}, B::Operator{<:PauliStringTS}; maxlength::Int=1000)
      295      Ls = qubitsize(A)
      296      Ps = periodicflags(A)
      263 -    N = Ls[1]                                                 
      297  
      298      p1s, c1s = A.strings, A.coeffs
      299      p2s, c2s = B.strings, B.coeffs
     ...
      303      length(p2s) == length(c2s) || throw(DimensionMismatch("str
           ings and coefficients must have the same length"))
      304  
      305      @inbounds for (p1, c1) in zip(p1s, c1s)
      273 -        rep1 = representative(p1)                             
      274 -        u1 = rep1.v | rep1.w                                  
      275 -        iszero(u1) && continue                                
      306 +        hrep = representative(p1)                             
      307          for (p2, c2) in zip(p2s, c2s)
      277 -            rep2 = representative(p2)                         
      278 -            u2 = rep2.v | rep2.w                              
      279 -            iszero(u2) && continue                            
      280 -                                                              
      281 -            seen = zero(typeof(u2))                           
      282 -            x1 = u1                                           
      283 -            while !iszero(x1)                                 
      284 -                i1 = trailing_zeros(x1) + 1                   
      285 -                x1 &= x1 - one(x1)                            
      286 -                                                              
      287 -                x2 = u2                                       
      288 -                while !iszero(x2)                             
      289 -                    i2 = trailing_zeros(x2) + 1               
      290 -                    x2 &= x2 - one(x2)                        
      291 -                                                              
      292 -                    s = mod(i1 - i2, N)                       
      293 -                    iszero(s) && (s = N)                      
      294 -                    bit = one(typeof(seen)) << (s - 1)        
      295 -                    (seen & bit != 0) && continue             
      296 -                    seen |= bit                               
      297 -                                                              
      298 -                    p, k = commutator(rep1, shift(rep2, Ls, Ps
          -, (s,)))                                                      
      299 -                    c = c1 * c2 * k                           
      300 -                    if (k != 0) && pauli_weight(p) < maxlength
      301 -                        setwith!(+, d, PauliStringTS{Ls,Ps}(p)
          -, c)                                                          
      302 -                    end                                       
      303 -                end                                           
      304 -            end                                               
      308 +            _orbit_edges_1d!(d, hrep, c1, p2, c2, Ls, Ps; maxl
          +ength=maxlength)                                              
      309          end
      310      end
      311      return d

  Ran 1 shell command (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/orbit_edge_stats.jl)
  ⎿  Wrote 71 lines to
     ../../../.claude/jobs/513a65b3/tmp/orbit_edge_stats.jl
       1 using PauliStrings
       2 using Dictionaries
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
       8 function Xtot(N)
       9     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
      10 end
     … +70 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● The edge-emitter refactor is correct but the timing regressed badly; that’s an
  implementation issue, not theory. I introduced abstract argument types in the
  hot helper (PauliString / PauliStringTS), which likely made the tight loop
  type-unstable. I’ll fix that immediately with parametric concrete signatures
  and re-measure.

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 1 line
      257      return op
      258  end
      259  
      260 -function _orbit_edges_1d!(d, hrep::PauliString, hcoeff, q::Pau
          -liStringTS, qcoeff, Ls::Tuple, Ps::Tuple; maxlength::Int=1000)
      260 +function _orbit_edges_1d!(d, hrep::P, hcoeff, q::Q, qcoeff, Ls
          +::Tuple, Ps::Tuple; maxlength::Int=1000) where {P<:PauliString
          +,Q<:PauliStringTS}                                            
      261      N = Ls[1]
      262      qrep = representative(q)
      263      uh = hrep.v | hrep.w

  Ran 1 shell command (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 1 line
      257      return op
      258  end
      259  
      260 -function _orbit_edges_1d!(d, hrep::P, hcoeff, q::Q, qcoeff, Ls
          -::Tuple, Ps::Tuple; maxlength::Int=1000) where {P<:PauliString
          -,Q<:PauliStringTS}                                            
      260 +function _orbit_edges_1d!(d, hrep::P, hcoeff, q::PauliStringTS
          +{Ls,Ps}, qcoeff; maxlength::Int=1000) where {P<:PauliString,Ls
          +,Ps}                                                          
      261      N = Ls[1]
      262      qrep = representative(q)
      263      uh = hrep.v | hrep.w

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 1 line
      305      @inbounds for (p1, c1) in zip(p1s, c1s)
      306          hrep = representative(p1)
      307          for (p2, c2) in zip(p2s, c2s)
      308 -            _orbit_edges_1d!(d, hrep, c1, p2, c2, Ls, Ps; maxl
          -ength=maxlength)                                              
      308 +            _orbit_edges_1d!(d, hrep, c1, p2, c2; maxlength=ma
          +xlength)                                                      
      309          end
      310      end
      311      return d

  Ran 1 shell command (ctrl+o to expand)

● Done. I refactored the orbit Liouvillian into an explicit per-orbit edge
  emitter backend and fixed the type-instability regression.

  Implemented

  Added internal edge emitter:

  _orbit_edges_1d!(d, hrep, hcoeff, q, qcoeff; maxlength)

  This emits the outgoing Pauli-orbit graph edges for one input TS Pauli orbit q
  under one Hamiltonian representative hrep.

  It is now used by:

  _commutator_1d_small_left!(...)
  orbit_liouvillian(...)

  So the backend is now structurally:

  for each input orbit q:
      emit graph edges q -> r
      accumulate into output dict

  rather than being just an opaque binary commutator loop.

  I also fixed the helper signature to carry Ls, Ps as type parameters:

  function _orbit_edges_1d!(d, hrep::P, hcoeff, q::PauliStringTS{Ls,Ps}, qcoeff;
  ...)

  The previous abstract signature caused a serious regression.

  Correctness

  Validated:

  orbit_liouvillian(Hts, Ots) == 1im * commutator(Hts, Ots)

  and against dense/resummed commutator for small systems:

  N=6 err_ref=0.0 err_full=0.0
  N=8 err_ref=0.0 err_full=0.0

  Edge statistics

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/orbit_edge_stats.jl

  N=32, M=12, lenO=4096

  full shifts over all orbits: 393216
  candidate active shifts:     76775  (19.52%)
  nonzero edges:               45569
  avg folded edges per input/orbit: 3.708
  max folded per input/orbit: 8

  N=64, M=14, lenO=16384

  full shifts over all orbits: 3145728
  candidate active shifts:     350724 (11.15%)
  nonzero edges:               208410
  avg folded edges per input/orbit: 4.240
  max folded per input/orbit: 9

  This is exactly the graph picture we wanted: the orbit graph degree is small,
  even though the full translation scan is large.

  Instrumentation after fix

  N=32, M=12

  ┌────────────────────────┬─────────┬─────────┬────────────────────────────┐
  │         Method         │  Time   │  Alloc  │         Main count         │
  ├────────────────────────┼─────────┼─────────┼────────────────────────────┤
  │ RK4                    │ 0.0167s │ 18.1 MB │           4 TS commutators │
  ├────────────────────────┼─────────┼─────────┼────────────────────────────┤
  │ Expanded Trotter       │ 0.0510s │ 57.2 MB │                  191 gates │
  ├────────────────────────┼─────────┼─────────┼────────────────────────────┤
  │ TrotterTS edge backend │ 0.0214s │ 28.8 MB │ 10 orbit Liouvillian calls │
  └────────────────────────┴─────────┴─────────┴────────────────────────────┘

  N=64, M=14

  ┌────────────────────────┬────────┬──────────┬────────────────────────────┐
  │         Method         │  Time  │  Alloc   │         Main count         │
  ├────────────────────────┼────────┼──────────┼────────────────────────────┤
  │ RK4                    │ 0.121s │  86.2 MB │           4 TS commutators │
  ├────────────────────────┼────────┼──────────┼────────────────────────────┤
  │ Expanded Trotter       │ 0.287s │ 470.7 MB │                  383 gates │
  ├────────────────────────┼────────┼──────────┼────────────────────────────┤
  │ TrotterTS edge backend │ 0.118s │ 138.3 MB │ 10 orbit Liouvillian calls │
  └────────────────────────┴────────┴──────────┴────────────────────────────┘

  Note: this instrumentation still labels the third row as Taylor2 because the
  script is old. In current code, public TrotterTS() is Krylov-only; the
  script’s local counted routine still uses a Taylor-style diagnostic to measure
  orbit-liouvillian cost. The edge backend itself is what matters here.

  Important outcome

  We now have:

  orbit graph degree ~4 per input/orbit

  for the tested evolved states.

  That means local-component exponentials are plausible: the graph is sparse
  enough that we should test closure growth under repeated orbit edges. If
  connected components remain modest under one orbit generator, exact local
  exponentiation could be a real TS-native primitive and avoid the expensive
  global Krylov orthogonalization.

✻ Sautéed for 3m 46s

❯ I see that you are using a 1d specializer. If I understand correctly this     
  will only work for the specific target. If this is the case, remove it from   
  the codebase immediately, to prevent further confusion                        

● You’re correct. The current active-shift backend is a 1D periodic specializer,
  so it risks optimizing for the target example rather than solving the general
  TS problem. I’ll remove it from the codebase now and leave only a generic
  orbit_liouvillian wrapper.

  Read 1 file (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Removed 70 lines
      257      return op
      258  end
      259  
      260 -function _orbit_edges_1d!(d, hrep::P, hcoeff, q::PauliStringTS
          -{Ls,Ps}, qcoeff; maxlength::Int=1000) where {P<:PauliString,Ls
          -,Ps}                                                          
      261 -    N = Ls[1]                                                 
      262 -    qrep = representative(q)                                  
      263 -    uh = hrep.v | hrep.w                                      
      264 -    uq = qrep.v | qrep.w                                      
      265 -    (iszero(uh) || iszero(uq)) && return d                    
      266 -                                                              
      267 -    seen = zero(typeof(uq))                                   
      268 -    xh = uh                                                   
      269 -    while !iszero(xh)                                         
      270 -        ih = trailing_zeros(xh) + 1                           
      271 -        xh &= xh - one(xh)                                    
      272 -                                                              
      273 -        xq = uq                                               
      274 -        while !iszero(xq)                                     
      275 -            iq = trailing_zeros(xq) + 1                       
      276 -            xq &= xq - one(xq)                                
      277 -                                                              
      278 -            s = mod(ih - iq, N)                               
      279 -            iszero(s) && (s = N)                              
      280 -            bit = one(typeof(seen)) << (s - 1)                
      281 -            (seen & bit != 0) && continue                     
      282 -            seen |= bit                                       
      283 -                                                              
      284 -            p, k = commutator(hrep, shift(qrep, Ls, Ps, (s,)))
      285 -            c = hcoeff * qcoeff * k                           
      286 -            if (k != 0) && pauli_weight(p) < maxlength        
      287 -                setwith!(+, d, PauliStringTS{Ls,Ps}(p), c)    
      288 -            end                                               
      289 -        end                                                   
      290 -    end                                                       
      291 -    return d                                                  
      292 -end                                                           
      293 -                                                              
      294 -function _commutator_1d_small_left!(d, A::Operator{<:PauliStri
          -ngTS}, B::Operator{<:PauliStringTS}; maxlength::Int=1000)     
      295 -    Ls = qubitsize(A)                                         
      296 -    Ps = periodicflags(A)                                     
      297 -                                                              
      298 -    p1s, c1s = A.strings, A.coeffs                            
      299 -    p2s, c2s = B.strings, B.coeffs                            
      300 -                                                              
      301 -    # check lengths to safely use `@inbounds`                 
      302 -    length(p1s) == length(c1s) || throw(DimensionMismatch("str
          -ings and coefficients must have the same length"))            
      303 -    length(p2s) == length(c2s) || throw(DimensionMismatch("str
          -ings and coefficients must have the same length"))            
      304 -                                                              
      305 -    @inbounds for (p1, c1) in zip(p1s, c1s)                   
      306 -        hrep = representative(p1)                             
      307 -        for (p2, c2) in zip(p2s, c2s)                         
      308 -            _orbit_edges_1d!(d, hrep, c1, p2, c2; maxlength=ma
          -xlength)                                                      
      309 -        end                                                   
      310 -    end                                                       
      311 -    return d                                                  
      312 -end                                                           
      313 -                                                              
      314 -function _binary_kernel_commutator_1d_small_left(A::Operator{<
          -:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real=0
          -, maxlength::Int=1000)                                        
      315 -    d = emptydict(A)                                          
      316 -    _commutator_1d_small_left!(d, A, B; maxlength=maxlength)  
      317 -    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      318 -    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilo
          -n) : o                                                        
      319 -end                                                           
      320 -                                                              
      260  function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
           gth::Int=1000)
      261      checklength(A, B)
      323 -    Ls = qubitsize(A)                                         
      324 -    Ps = periodicflags(A)                                     
      325 -    if length(Ls) == 1 && Ps[1] && length(A) <= 8 && qubitleng
          -th(A) <= 8*sizeof(representative(A.strings[1]).v)             
      326 -        d = emptydict(A)                                      
      327 -        _commutator_1d_small_left!(d, A, B; maxlength=maxlengt
          -h)                                                            
      328 -        o = typeof(A)(collect(keys(d)), collect(values(d)))   
      329 -        o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsil
          -on) : o                                                       
      330 -        return 1im * o / hbar                                 
      331 -    end                                                       
      262      return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
           axlength) / hbar
      263  end
      264  

● Update(src/translation_symmetry.jl)
  ⎿  Removed 3 lines
      266      checklength(A, B)
      267      Ls = qubitsize(A)
      268      Ps = periodicflags(A)
      269 -    if op === commutator && length(Ls) == 1 && Ps[1] && length
          -(A) <= 8 && qubitlength(A) <= 8*sizeof(representative(A.string
          -s[1]).v)                                                      
      270 -        return _binary_kernel_commutator_1d_small_left(A, B; e
          -psilon=epsilon, maxlength=maxlength)                          
      271 -    end                                                       
      269  
      270      d = emptydict(A)
      271      p1s, c1s = A.strings, A.coeffs

  Searched for 1 pattern, ran 1 shell command (ctrl+o to expand)

● Removed the 1D specializer from the codebase.

  Removed

  Deleted all of:

  _orbit_edges_1d!
  _commutator_1d_small_left!
  _binary_kernel_commutator_1d_small_left

  and removed the 1D dispatch branch from binary_kernel.

  Also simplified orbit_liouvillian back to the generic implementation:

  function orbit_liouvillian(A::Operator{<:PauliStringTS},
  B::Operator{<:PauliStringTS}; hbar=1, epsilon=0, maxlength=1000)
      checklength(A, B)
      return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
  end

  Verified

  Searched for leftover specializer names:

  1d
  small_left
  orbit_edges
  _binary_kernel_commutator

  No matches remain in src/translation_symmetry.jl or src/evolution.jl.

  test/evolution.jl passes:

  evolve | 8 passed

  So the codebase no longer contains the target-specific 1D active-shift
  optimization.

✻ Cooked for 1m 25s

❯ now check performance                                                         

  Ran 2 shell commands (ctrl+o to expand)

● Performance after removing the 1D specializer and making TrotterTS()
  Krylov-only is bad, as expected.

  Benchmark: MFIM, dt=0.1, no fout.

  ┌─────────────────┬────────────┬──────────────┬──────────────────────────┐
  │      Case       │    RK4     │  Expanded    │         Generic          │
  │                 │            │   Trotter    │   TrotterTS(:krylov,4)   │
  ├─────────────────┼────────────┼──────────────┼──────────────────────────┤
  │ N=32, M=12,     │    1.19s,  │ 1.28s, 2.38  │           5.23s, 9.68 GB │
  │ steps=50        │    1.07 GB │           GB │                          │
  ├─────────────────┼────────────┼──────────────┼──────────────────────────┤
  │ N=64, M=14,     │    2.27s,  │ 4.01s, 8.07  │           5.88s, 8.85 GB │
  │ steps=20        │    1.35 GB │           GB │                          │
  └─────────────────┴────────────┴──────────────┴──────────────────────────┘

  Interpretation

  Generic Krylov-orbit TrotterTS is currently slower than both RK4 and expanded
  Trotter.

  This is not surprising because:

  TrotterTS(:krylov,4)

  does many TS Liouvillian applications plus orthogonalization and operator
  combinations per orbit factor.

  For MFIM:

  orbit factors per step = 5
  krylovdim = 4
  Liouvillian applications ≈ 20 per step

  RK4 does:

  4 full-H commutators per step

  Each full-H commutator has 3 orbit terms, so roughly:

  ~12 orbit-term commutator equivalents per step

  So the current generic Krylov implementation does more operator work, plus
  more allocation-heavy orthogonalization.

  What this tells us

  Removing the 1D specializer makes the state honest:

  - We have a conceptually proper TS orbit-Trotter method.
  - But the generic backend is not performant.
  - The previous speedups came partly from the 1D active-shift specializer and
  partly from Taylor-2, both of which we intentionally removed/rejected.

  The real next milestone remains:

  ▎ a generic dimension-independent orbit graph edge generator / local-component
  ▎ exponential backend.

  Not a 1D hot path. Not Taylor-2.

  Without that, Krylov is too expensive.

✻ Crunched for 2m 4s

❯ /memory 
  ⎿  Opened memory file at ./CLAUDE.md

     ▎ Using $EDITOR="nvim". To change editor, set $EDITOR or $VISUAL environment variable.
❯ /memory 
  ⎿  Opened memory file at ~/.claude/CLAUDE.md

     ▎ Using $EDITOR="nvim". To change editor, set $EDITOR or $VISUAL environment variable.

❯ /memory 
  ⎿  Cancelled memory editing

❯ go                                                                            

● I’ll continue by adding a generic, dimension-independent orbit-graph edge
  interface, but I won’t route the hot evolution through it until it is
  measured. That keeps the code honest: the graph primitive exists for
  validation/local-component experiments without sneaking in target-specific
  performance paths.

  Read 3 files, ran 1 shell command (ctrl+o to expand)

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line, removed 1 line
      23  export k_local_basis, z_basis, x_basis, y_basis, complete_basis
      24  export Circuits
      25  export pauli_weight, support, pauli_rotation
      26 -export TrotterGate, trotterize, trotter_step!, orbit_liouvillia
         -n                                                              
      26 +export TrotterGate, trotterize, trotter_step!, orbit_edges, orb
         +it_liouvillian                                                 
      27  export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
          r, TrotterTS, RK4, DOPRI5, Exact
      28  
      29  using Random

● Update(src/translation_symmetry.jl)
  ⎿  Added 29 lines, removed 1 line
      257      return op
      258  end
      259  
      260 +function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStr
          +ingTS; epsilon::Real=0, maxlength::Int=1000)                  
      261 +    checklength(A, q)                                         
      262 +    Ls = qubitsize(A)                                         
      263 +    Ps = periodicflags(A)                                     
      264 +    d = emptydict(A)                                          
      265 +    qrep = representative(q)                                  
      266 +    for (p, c) in zip(A.strings, A.coeffs)                    
      267 +        prep = representative(p)                              
      268 +        for s in all_shifts(paulistringtype(A))               
      269 +            out, k = commutator(prep, shift(qrep, Ls, Ps, s)) 
      270 +            coeff = c * k                                     
      271 +            if (k != 0) && (abs(coeff) > epsilon) && pauli_wei
          +ght(out) < maxlength                                          
      272 +                setwith!(+, d, PauliStringTS{Ls,Ps}(out), coef
          +f)                                                            
      273 +            end                                               
      274 +        end                                                   
      275 +    end                                                       
      276 +    return typeof(A)(collect(keys(d)), collect(values(d)))    
      277 +end                                                           
      278 +                                                              
      279  function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
           gth::Int=1000)
      280      checklength(A, B)
      262 -    return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
          -axlength) / hbar                                              
      281 +    d = emptydict(A)                                          
      282 +    for (q, cq) in zip(B.strings, B.coeffs)                   
      283 +        edges = orbit_edges(A, q; epsilon=epsilon, maxlength=m
          +axlength)                                                     
      284 +        for (p, c) in zip(edges.strings, edges.coeffs)        
      285 +            setwith!(+, d, p, cq * c)                         
      286 +        end                                                   
      287 +    end                                                       
      288 +    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      289 +    o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) 
          +: o                                                           
      290 +    return 1im * o / hbar                                     
      291  end
      292  
      293  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)

  Ran 1 shell command (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/generic_orbit_edge_stats.jl)
  ⎿  Wrote 35 lines to
     ../../../.claude/jobs/513a65b3/tmp/generic_orbit_edge_stats.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +29 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Implemented a generic dimension-independent edge interface.

  Added

  Exported:

  orbit_edges(A::OperatorTS, q::PauliStringTS; epsilon=0, maxlength=1000)

  This returns the outgoing Pauli-orbit graph edges for one input orbit q under
  orbit/operator A.

  Mathematically:

  q \mapsto \{(r, c_{r,q})\}

  where:

  L_A([q]) = i[A, [q]] = i \sum_r c_{r,q} [r].

  Implementation is generic: it uses existing all_shifts, shift, commutator, and
  TS canonicalization. No 1D specializer, no target-specific path.

  Also refactored:

  orbit_liouvillian(A, B)

  to explicitly sum orbit_edges(A, q) over the strings in B.

  So the code now exposes the graph structure, albeit through a straightforward
  generic implementation.

  Correctness

  Validated:

  orbit_liouvillian(Hts, Ots) == 1im * commutator(Hts, Ots)

  and for a single input orbit:

  orbit_edges(Hts, q) == commutator(Hts, OperatorTS(q))

  up to zero resummed norm error for small systems:

  N=6 err=0.0 edge_err=0.0
  N=8 err=0.0 edge_err=0.0

  Generic graph statistics

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/generic_orbit_edge_stats.jl

  N=32, M=12, lenO=4096

  ALL folded_edges = 45569
  avg per input/orbit = 3.708
  max = 8

  Per orbit:

  orbit 1 avg=3.958 max=8
  orbit 2 avg=3.497 max=6
  orbit 3 avg=3.670 max=7

  N=64, M=14, lenO=16384

  ALL folded_edges = 208410
  avg per input/orbit = 4.240
  max = 9

  Per orbit:

  orbit 1 avg=3.972 max=8
  orbit 2 avg=4.157 max=8
  orbit 3 avg=4.592 max=9

  Why this matters

  The graph degree is small:

  ~4 folded outgoing edges per input orbit

  That is the theoretical signal we wanted.

  Even though generic orbit_edges still scans all translations internally, the
  folded orbit graph itself is sparse. This makes local-component exponentials
  or a better graph mul! plausible.

  Caveat

  The current generic orbit_edges allocates and scans all shifts. It is not
  intended as the final fast backend.

  Its purpose is to make the graph explicit and measurable without 1D-specific
  optimization. Next step should be:

  build local connected components using orbit_edges
  measure closure growth

  If components remain small, local-component exponentials are the next serious
  candidate for a true TS-native orbit flow.

✻ Worked for 2m 9s

❯ continue                                                                      

● I’ll measure local orbit-graph closure growth next. That tells us whether
  exact local-component exponentials are plausible before implementing them.

● Write(~/.claude/jobs/513a65b3/tmp/orbit_component_growth.jl)
  ⎿  Wrote 76 lines to
     ../../../.claude/jobs/513a65b3/tmp/orbit_component_growth.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +76 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I measured orbit-graph component growth using the generic orbit_edges
  interface.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/orbit_component_growth.jl

  This samples TS strings from evolved states and repeatedly closes them under
  one Hamiltonian orbit’s edge map:

  seed q
  depth 1: orbit_edges(Ha, q)
  depth 2: orbit_edges(Ha, all new nodes)
  ...

  No 1D-specialized code involved in the implementation; this uses generic
  orbit_edges.

  Results

  N=16, M=10, lenO=1024, depth 5

  ┌───────┬────────────────────┬─────┬─────┬─────┐
  │ Orbit │ Avg component size │ p50 │ p90 │ Max │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     1 │               10.0 │   4 │  16 │  16 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     2 │               12.6 │   8 │  16 │  63 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     3 │               13.7 │   8 │  32 │  32 │
  └───────┴────────────────────┴─────┴─────┴─────┘

  N=32, M=12, lenO=4096, depth 4

  ┌───────┬────────────────────┬─────┬─────┬─────┐
  │ Orbit │ Avg component size │ p50 │ p90 │ Max │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     1 │               13.6 │  16 │  31 │  31 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     2 │               14.2 │  16 │  16 │  57 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     3 │               14.4 │   8 │  31 │  31 │
  └───────┴────────────────────┴─────┴─────┴─────┘

  N=64, M=14, lenO=16384, depth 3

  ┌───────┬────────────────────┬─────┬─────┬─────┐
  │ Orbit │ Avg component size │ p50 │ p90 │ Max │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     1 │               15.6 │  15 │  42 │  42 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     2 │                8.4 │   8 │  15 │  15 │
  ├───────┼────────────────────┼─────┼─────┼─────┤
  │     3 │               22.6 │  15 │  64 │  64 │
  └───────┴────────────────────┴─────┴─────┴─────┘

  Interpretation

  This is encouraging.

  The orbit graph does not immediately explode for sampled states. Components
  are often tiny:

  ~4, 8, 15, 16, 31, 64 nodes

  Those are sizes where local dense exponentiation of a component is plausible.

  This suggests a possible true TS-native orbit exponential primitive:

  for each orbit flow e^{τL_a}:
      decompose active TS Pauli graph into connected components under L_a
      exponentiate each component locally
      update coefficients

  This would be conceptually different from RK4:

  RK4: polynomial in full-H Liouvillian
  Orbit Trotter: product of local exact exponentials of orbit Liouvillian
  components

  and it would be a genuine Pauli-orbit propagation rule.

  Caveats

  - These are sampled components, not exhaustive over all active strings.
  - The closure depths were finite:
    - depth 5 for N=16,
    - depth 4 for N=32,
    - depth 3 for N=64.
  - Some components may keep growing at greater depth.
  - The generic orbit_edges currently scans all shifts and allocates, so
  component discovery itself is not yet optimized.

  Next milestone

  Implement a prototype for one orbit flow:

  orbit_flow_component_exp(Ha, O, τ)

  Algorithm:

  1. For each active string in O, build its connected component under
  orbit_edges(Ha, q).
  2. Build the finite Liouvillian matrix on that component.
  3. Apply:

  exp(τ * A_component)

  3. to coefficients in that component.
  4. Reassemble OperatorTS.

  Then validate against dense exact for one orbit flow at small N:

  e^{iτ Ha} O e^{-iτ Ha}

  If component closure remains small and exact component exponentials are not
  too expensive, this is the first real candidate for a TS-native Trotter
  primitive.

✻ Baked for 1m 1s

❯ continue                                                                      

● Update(src/evolution.jl)
  ⎿  Added 9 lines, removed 5 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, substep=:krylov, krylovdim=4)         
      29 +    TrotterTS(; order=2, substep=:krylov, krylovdim=4, componen
         +tmax=256)                                                      
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34 -Liouvillians. Each orbit flow is applied as a matrix-free Arnol
         -di/Krylov exponential                                          
      35 -action with dimension `krylovdim`.                             
      34 +Liouvillians. `substep=:krylov` applies each orbit flow as a ma
         +trix-free                                                      
      35 +Arnoldi/Krylov exponential action with dimension `krylovdim`. `
         +substep=:component`                                            
      36 +builds connected components of the Pauli-orbit graph up to `com
         +ponentmax` nodes and                                           
      37 +applies a dense local exponential on each component, falling ba
         +ck to Krylov when a                                            
      38 +component exceeds the cap.                                     
      39  """
      40  struct TrotterTS <: AbstractEvolutionMethod
      41      order::Int
      42      substep::Symbol
      43      krylovdim::Int
      44 +    componentmax::Int                                          
      45  end
      42 -TrotterTS(; order::Integer=2, substep::Symbol=:krylov, krylovdi
         -m::Integer=4) =                                                
      43 -    TrotterTS(Int(order), substep, Int(krylovdim))             
      46 +TrotterTS(; order::Integer=2, substep::Symbol=:krylov, krylovdi
         +m::Integer=4, componentmax::Integer=256) =                     
      47 +    TrotterTS(Int(order), substep, Int(krylovdim), Int(componen
         +tmax))                                                         
      48  
      49  """
      50      RK4()

● I’ll add a prototype substep=:component orbit flow: build connected components
  with generic orbit_edges, exponentiate each component densely, and fall back
  to Krylov if the component exceeds a cap. Then I’ll validate on small systems
  before benchmarking.

● Update(src/evolution.jl)
  ⎿  Added 70 lines, removed 3 lines
      342      return truncation(_axpy_combination(V, coeffs, O))
      343  end
      344  
      345 +function _build_orbit_component(Ha::Operator{<:PauliStringTS},
          + seed::PauliStringTS, componentmax::Int, maxlength::Int)      
      346 +    component = typeof(seed)[seed]                            
      347 +    index = Dict{typeof(seed),Int}(seed => 1)                 
      348 +    queue_index = 1                                           
      349 +    capped = false                                            
      350 +    while queue_index <= length(component)                    
      351 +        q = component[queue_index]                            
      352 +        queue_index += 1                                      
      353 +        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      354 +        for r in edges.strings                                
      355 +            if !haskey(index, r)                              
      356 +                if length(component) >= componentmax          
      357 +                    capped = true                             
      358 +                    return component, index, capped           
      359 +                end                                           
      360 +                push!(component, r)                           
      361 +                index[r] = length(component)                  
      362 +            end                                               
      363 +        end                                                   
      364 +    end                                                       
      365 +    return component, index, capped                           
      366 +end                                                           
      367 +                                                              
      368 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, hbar::Real, componentm
          +ax::Int, krylovdim::Int, truncation)                          
      369 +    componentmax > 0 || throw(ArgumentError("componentmax must
          + be positive"))                                               
      370 +    maxlength = 1000                                          
      371 +    out = zero(O)                                             
      372 +    assigned = Set{eltype(O.strings)}()                       
      373 +                                                              
      374 +    for q0 in O.strings                                       
      375 +        (q0 in assigned) && continue                          
      376 +        component, index, capped = _build_orbit_component(Ha, 
          +q0, componentmax, maxlength)                                  
      377 +        if capped                                             
      378 +            seed_op = typeof(O)([q0], [O.coeffs[findfirst(==(q
          +0), O.strings)]])                                             
      379 +            out = out + _orbit_flow_krylov(Ha, seed_op, dt, hb
          +ar, krylovdim, truncation)                                    
      380 +            push!(assigned, q0)                               
      381 +            continue                                          
      382 +        end                                                   
      383 +                                                              
      384 +        ncomp = length(component)                             
      385 +        A = zeros(ComplexF64, ncomp, ncomp)                   
      386 +        for (j, q) in enumerate(component)                    
      387 +            edges = orbit_edges(Ha, q; hbar=hbar, maxlength=ma
          +xlength)                                                      
      388 +            for (r, c) in zip(edges.strings, edges.coeffs)    
      389 +                i = get(index, r, 0)                          
      390 +                i != 0 && (A[i, j] += 1im * c / hbar)         
      391 +            end                                               
      392 +        end                                                   
      393 +                                                              
      394 +        coeffs = zeros(ComplexF64, ncomp)                     
      395 +        for (q, c) in zip(O.strings, O.coeffs)                
      396 +            j = get(index, q, 0)                              
      397 +            if j != 0                                         
      398 +                coeffs[j] += c                                
      399 +                push!(assigned, q)                            
      400 +            end                                               
      401 +        end                                                   
      402 +        coeffs2 = exp(dt * A) * coeffs                        
      403 +        out = out + typeof(O)(component, coeffs2)             
      404 +    end                                                       
      405 +    return truncation(out)                                    
      406 +end                                                           
      407 +                                                              
      408  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
           runcation)
      346 -    method.substep === :krylov ||                             
      347 -        throw(ArgumentError("unknown TrotterTS substep $(metho
          -d.substep); only :krylov is currently implemented"))          
      348 -    return _orbit_flow_krylov(Ha, O, dt, hbar, method.krylovdi
          -m, truncation)                                                
      409 +    if method.substep === :krylov                             
      410 +        return _orbit_flow_krylov(Ha, O, dt, hbar, method.kryl
          +ovdim, truncation)                                            
      411 +    elseif method.substep === :component                      
      412 +        return _orbit_flow_component(Ha, O, dt, hbar, method.c
          +omponentmax, method.krylovdim, truncation)                    
      413 +    else                                                      
      414 +        throw(ArgumentError("unknown TrotterTS substep $(metho
          +d.substep); expected :krylov or :component"))                 
      415 +    end                                                       
      416  end
      417  
      418  function _orbit_terms(H::Operator{<:PauliStringTS})

  Ran 1 shell command (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      384          ncomp = length(component)
      385          A = zeros(ComplexF64, ncomp, ncomp)
      386          for (j, q) in enumerate(component)
      387 -            edges = orbit_edges(Ha, q; hbar=hbar, maxlength=ma
          -xlength)                                                      
      387 +            edges = orbit_edges(Ha, q; maxlength=maxlength)   
      388              for (r, c) in zip(edges.strings, edges.coeffs)
      389                  i = get(index, r, 0)
      390                  i != 0 && (A[i, j] += 1im * c / hbar)

  Ran 1 shell command (ctrl+o to expand)

● Implemented a prototype local-component orbit exponential substep.

  Added

  TrotterTS now supports:

  TrotterTS(substep=:component, componentmax=256)

  alongside:

  TrotterTS(substep=:krylov, krylovdim=4)

  Component substep algorithm

  For each orbit flow:

  e^{τ L_a} O

  it now can:

  1. Build connected components of the Pauli-orbit graph using generic:

  orbit_edges(Ha, q)
  2. For each component, build a dense matrix representation of the orbit
  Liouvillian:

  A_{ij}
  3. Apply exact local dense exponential:

  exp(τ * A_component) * coeffs
  4. Reassemble OperatorTS.

  If a component exceeds componentmax, it falls back to Krylov for that seed.

  This is still generic and dimension-independent. No 1D specializer.

  Correctness smoke test

  For N=6, TFIM, dt=0.05, compared trajectory against dense exact:

  TrotterTS(substep=:component, componentmax=128):
    err = 2.766e-4

  TrotterTS(substep=:krylov, krylovdim=4):
    err = 2.766e-4

  The two agree to displayed precision, which is expected if components are
  closed and exponentiated exactly.

  This is a good sign: component exponentials are reproducing the intended
  orbit-flow exponential, not Taylor/RK-like behavior.

  Caveat

  This is a first prototype and likely allocates heavily:

  - uses Dict,
  - builds components repeatedly,
  - dense exp per component,
  - not yet caching components or matrices.

  So the next step is to benchmark and instrument it, especially:

  number of components
  component sizes
  how often componentmax fallback happens
  time in component building vs matrix exp
  allocations

  But conceptually, this is now the correct TS-native Trotter primitive
  direction:

  Strang over orbit Hamiltonians
  +
  exact local exponentials on Pauli-orbit graph components

✻ Cooked for 2m 45s

❯ ok go ahead                                                                   
● Write(~/.claude/jobs/513a65b3/tmp/component_perf_stats.jl)
  ⎿  Wrote 74 lines to
     ../../../.claude/jobs/513a65b3/tmp/component_perf_stats.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +80 lines (ctrl+o to expand)

● I’ll instrument component flow now: component sizes/fallbacks, allocations,
  and timings versus Krylov/RK4/expanded Trotter on smaller cases first so we
  don’t waste time if the prototype is too allocation-heavy.

  Ran 1 shell command (ctrl+o to expand)

● I benchmarked and instrumented the local-component exponential prototype.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/component_perf_stats.jl

  Component stats

  N=16, M=10, steps=20

  Representative state after grow:

  lenO = 1024
  componentmax = 128

  ┌───────┬────────────┬────────┬───────────────────────┬─────┬─────┬─────┐
  │ Orbit │ Components │ Capped │ Total component nodes │ p50 │ p90 │ Max │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     1 │        229 │      0 │                  2307 │   4 │  32 │ 128 │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     2 │        217 │      0 │                  1555 │   4 │  16 │  64 │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     3 │        251 │      0 │                  2181 │   8 │  16 │  64 │
  └───────┴────────────┴────────┴───────────────────────┴─────┴─────┴─────┘

  N=32, M=12, steps=10

  Representative state after grow:

  lenO = 4096
  componentmax = 128

  ┌───────┬────────────┬────────┬───────────────────────┬─────┬─────┬─────┐
  │ Orbit │ Components │ Capped │ Total component nodes │ p50 │ p90 │ Max │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     1 │        782 │      0 │                  9216 │   8 │  32 │ 128 │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     2 │        697 │      0 │                  7003 │   4 │  16 │  64 │
  ├───────┼────────────┼────────┼───────────────────────┼─────┼─────┼─────┤
  │     3 │        670 │      3 │                  9911 │   8 │  32 │ 128 │
  └───────┴────────────┴────────┴───────────────────────┴─────┴─────┴─────┘

  So the graph component idea is structurally real:

  - components are mostly small,
  - many p50/p90 sizes are reasonable,
  - few caps at N=32.

  Performance

  N=16, M=10, steps=20

  ┌───────────────────────┬────────┬──────────┬───────────┐
  │        Method         │  Time  │  Alloc   │ Final len │
  ├───────────────────────┼────────┼──────────┼───────────┤
  │ RK4                   │ 0.057s │ 0.063 GB │      1024 │
  ├───────────────────────┼────────┼──────────┼───────────┤
  │ Expanded Trotter      │ 0.066s │ 0.116 GB │       514 │
  ├───────────────────────┼────────┼──────────┼───────────┤
  │ TrotterTS(:krylov,4)  │ 0.476s │ 1.082 GB │      1024 │
  ├───────────────────────┼────────┼──────────┼───────────┤
  │ TrotterTS(:component) │ 2.293s │ 4.811 GB │      1024 │
  └───────────────────────┴────────┴──────────┴───────────┘

  N=32, M=12, steps=10

  ┌───────────────────────┬────────┬───────────┬───────────┐
  │        Method         │  Time  │   Alloc   │ Final len │
  ├───────────────────────┼────────┼───────────┼───────────┤
  │ RK4                   │ 0.118s │  0.102 GB │      4096 │
  ├───────────────────────┼────────┼───────────┼───────────┤
  │ Expanded Trotter      │ 0.248s │  0.412 GB │      1401 │
  ├───────────────────────┼────────┼───────────┼───────────┤
  │ TrotterTS(:krylov,4)  │ 0.815s │  1.699 GB │      4096 │
  ├───────────────────────┼────────┼───────────┼───────────┤
  │ TrotterTS(:component) │ 9.004s │ 19.777 GB │      4096 │
  └───────────────────────┴────────┴───────────┴───────────┘

  Conclusion

  The theory is promising, but the naive component implementation is unusably
  slow.

  Why:

  - It rebuilds components every orbit flow and every time step.
  - It calls orbit_edges repeatedly.
  - It creates hundreds of dense matrices and computes hundreds of dense
  exponentials.
  - It allocates huge amounts of temporary data.

  This does not invalidate the Pauli-orbit component idea. It shows the naive
  implementation is the wrong one.

  What we learned

  Component closure is small enough to be interesting:

  p50 ~ 4–8
  p90 ~ 16–32

  But component exponentiation must be cached.

  The correct architecture is likely:

  precompute orbit graph components and local matrices for each Hamiltonian
  orbit
  reuse them across all time steps

  Since the orbit graph under a fixed Ha is independent of time, rebuilding it
  each step is wasteful.

  Next practical direction

  For a fixed truncation/candidate active basis, build an OrbitFlowPlan:

  struct OrbitFlowPlan
      components
      matrices
      expmats_for_dt
      index_maps
  end

  Then each orbit flow only does coefficient gather → matrix multiply → scatter.

  The open issue is basis adaptivity: trim changes active strings. But we can
  start with a dynamic plan that expands/caches components as new strings
  appear, rather than rebuilding from scratch.

  That would test whether the component-exponential idea can become practical.

✻ Cooked for 1m 35s

❯ Ok it seems you are confident with the theoretical advantage of our method.   
  before optimising the engineering problems, clean up the code to exactly and  
  concisely implement this method. fan out an agent to write the theory down    
  into a markdown file. be clear but not unnecessarily verbose. This will act   
  as the mathematical proposal document.                                        

● Agent(Write TS theory doc) Haiku 4.5
  ⎿  Backgrounded agent (↓ to manage · ctrl+o to expand)

● Agent "Write TS theory doc" completed · 0s

● I’ll clean the code down to the concise theoretical method: TrotterTS as
  Strang over orbit Liouvillian exponentials, with Krylov as the only subflow
  and orbit_edges/orbit_liouvillian as the generic graph primitives. I’ve also
  started a background agent to write the mathematical proposal markdown.

  Read 2 files (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 5 lines, removed 10 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, substep=:krylov, krylovdim=4, componen
         -tmax=256)                                                      
      29 +    TrotterTS(; order=2, krylovdim=4)                          
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34 -Liouvillians. `substep=:krylov` applies each orbit flow as a ma
         -trix-free                                                      
      35 -Arnoldi/Krylov exponential action with dimension `krylovdim`. `
         -substep=:component`                                            
      36 -builds connected components of the Pauli-orbit graph up to `com
         -ponentmax` nodes and                                           
      37 -applies a dense local exponential on each component, falling ba
         -ck to Krylov when a                                            
      38 -component exceeds the cap.                                     
      34 +Liouvillians. Each orbit flow is applied as a matrix-free Arnol
         +di/Krylov exponential                                          
      35 +action with dimension `krylovdim`.                             
      36  """
      37  struct TrotterTS <: AbstractEvolutionMethod
      38      order::Int
      42 -    substep::Symbol                                            
      39      krylovdim::Int
      44 -    componentmax::Int                                          
      40  end
      46 -TrotterTS(; order::Integer=2, substep::Symbol=:krylov, krylovdi
         -m::Integer=4, componentmax::Integer=256) =                     
      47 -    TrotterTS(Int(order), substep, Int(krylovdim), Int(componen
         -tmax))                                                         
      41 +TrotterTS(; order::Integer=2, krylovdim::Integer=4) =          
      42 +    TrotterTS(Int(order), Int(krylovdim))                      
      43  
      44  """
      45      RK4()

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 70 lines
      337      return truncation(_axpy_combination(V, coeffs, O))
      338  end
      339  
      340 -function _build_orbit_component(Ha::Operator{<:PauliStringTS},
          - seed::PauliStringTS, componentmax::Int, maxlength::Int)      
      341 -    component = typeof(seed)[seed]                            
      342 -    index = Dict{typeof(seed),Int}(seed => 1)                 
      343 -    queue_index = 1                                           
      344 -    capped = false                                            
      345 -    while queue_index <= length(component)                    
      346 -        q = component[queue_index]                            
      347 -        queue_index += 1                                      
      348 -        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      349 -        for r in edges.strings                                
      350 -            if !haskey(index, r)                              
      351 -                if length(component) >= componentmax          
      352 -                    capped = true                             
      353 -                    return component, index, capped           
      354 -                end                                           
      355 -                push!(component, r)                           
      356 -                index[r] = length(component)                  
      357 -            end                                               
      358 -        end                                                   
      359 -    end                                                       
      360 -    return component, index, capped                           
      361 -end                                                           
      362 -                                                              
      363 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, hbar::Real, componentm
          -ax::Int, krylovdim::Int, truncation)                          
      364 -    componentmax > 0 || throw(ArgumentError("componentmax must
          - be positive"))                                               
      365 -    maxlength = 1000                                          
      366 -    out = zero(O)                                             
      367 -    assigned = Set{eltype(O.strings)}()                       
      368 -                                                              
      369 -    for q0 in O.strings                                       
      370 -        (q0 in assigned) && continue                          
      371 -        component, index, capped = _build_orbit_component(Ha, 
          -q0, componentmax, maxlength)                                  
      372 -        if capped                                             
      373 -            seed_op = typeof(O)([q0], [O.coeffs[findfirst(==(q
          -0), O.strings)]])                                             
      374 -            out = out + _orbit_flow_krylov(Ha, seed_op, dt, hb
          -ar, krylovdim, truncation)                                    
      375 -            push!(assigned, q0)                               
      376 -            continue                                          
      377 -        end                                                   
      378 -                                                              
      379 -        ncomp = length(component)                             
      380 -        A = zeros(ComplexF64, ncomp, ncomp)                   
      381 -        for (j, q) in enumerate(component)                    
      382 -            edges = orbit_edges(Ha, q; maxlength=maxlength)   
      383 -            for (r, c) in zip(edges.strings, edges.coeffs)    
      384 -                i = get(index, r, 0)                          
      385 -                i != 0 && (A[i, j] += 1im * c / hbar)         
      386 -            end                                               
      387 -        end                                                   
      388 -                                                              
      389 -        coeffs = zeros(ComplexF64, ncomp)                     
      390 -        for (q, c) in zip(O.strings, O.coeffs)                
      391 -            j = get(index, q, 0)                              
      392 -            if j != 0                                         
      393 -                coeffs[j] += c                                
      394 -                push!(assigned, q)                            
      395 -            end                                               
      396 -        end                                                   
      397 -        coeffs2 = exp(dt * A) * coeffs                        
      398 -        out = out + typeof(O)(component, coeffs2)             
      399 -    end                                                       
      400 -    return truncation(out)                                    
      401 -end                                                           
      402 -                                                              
      340  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
           runcation)
      404 -    if method.substep === :krylov                             
      405 -        return _orbit_flow_krylov(Ha, O, dt, hbar, method.kryl
          -ovdim, truncation)                                            
      406 -    elseif method.substep === :component                      
      407 -        return _orbit_flow_component(Ha, O, dt, hbar, method.c
          -omponentmax, method.krylovdim, truncation)                    
      408 -    else                                                      
      409 -        throw(ArgumentError("unknown TrotterTS substep $(metho
          -d.substep); expected :krylov or :component"))                 
      410 -    end                                                       
      341 +    return _orbit_flow_krylov(Ha, O, dt, hbar, method.krylovdi
          +m, truncation)                                                
      342  end
      343  
      344  function _orbit_terms(H::Operator{<:PauliStringTS})

●

● Task Output(non-blocking) abada3702ccbbbeb9
  ⎿  Error: No task found with ID: abada3702ccbbbeb9

  Searched for 2 patterns, read 1 file (ctrl+o to expand)

● Write(benchmark/ts_native_trotter_theory.md)
  ⎿  Wrote 145 lines to benchmark/ts_native_trotter_theory.md
      1 # TS-native Trotter proposal
      2 
      3 ## Motivation
      4 
      5 The current translation-symmetric (TS) Trotter path expands
      6 
      7 ```math
      8 H = \sum_a c_a \sum_r T_r P_a T_r^\dagger
      9 ```
     … +168 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Cleaned up the implementation and wrote the theory proposal.

  Code cleanup

  Removed the experimental local-component substep from TrotterTS.

  TrotterTS is now concise again:

  TrotterTS(; order=2, krylovdim=4)

  with fields:

  order::Int
  krylovdim::Int

  It implements only:

  Strang over TS orbit Liouvillians
  +
  Krylov exponential action for each orbit flow

  No substep, no componentmax, no Taylor-2, no local-component prototype in
  production code.

  Kept the generic graph primitives:

  orbit_edges(A::OperatorTS, q::PauliStringTS)
  orbit_liouvillian(A::OperatorTS, B::OperatorTS)

  These are generic, dimension-independent, and not optimized for the target
  example.

  Theory document

  Created:

  benchmark/ts_native_trotter_theory.md

  It covers:

  - motivation,
  - orbit-level Strang splitting,
  - Pauli-orbit propagation rule,
  - how orbit_edges represents the graph,
  - difference from RK4 and expanded Trotter,
  - implementation strategy,
  - validation goals,
  - references.

  Tests

  Evolution tests pass:

  evolve | 8 passed

✻ Cooked for 2m 58s

※ recap: We’re designing a generic TS-native second-order Trotter method using orbit Liouvillian exponentials instead of expanded translated gates. Next, validate the cleaned Krylov-based implementation and theory doc, then decide whether to optimize orbit_edges.

❯ very good. now write a small paragraph about the other attempts we made with  
  taylor2 etc. This is only for tracking what we did                            

  Read 1 file (ctrl+o to expand)

● Update(benchmark/ts_native_trotter_theory.md)
  ⎿  Added 4 lines
      114  4. Investigate local connected components of the orbit graph. 
           If components are small, build local matrices and apply dense 
           exponentials per component.
      115  5. Avoid target-specific assumptions such as 1D nearest-neighb
           or structure or commuting MFIM orbits in the core method.
      116  
      117 +## Attempts tried during exploration                          
      118 +                                                              
      119 +We tried a cheap orbit-Strang variant where each orbit flow wa
          +s replaced by a second-order Taylor approximation, `O + τ L_a 
          +O + τ² L_a² O / 2`. This was useful diagnostically and showed 
          +that orbit-level splitting can reduce the expanded-gate work, 
          +but it is too close in spirit to a commutator-based integrator
          + and should not be considered the final TS Trotter primitive. 
          +We also prototyped local connected-component exponentials on t
          +he Pauli-orbit graph; small sampled components suggested theor
          +etical promise, but the naive implementation rebuilt component
          +s and dense matrices every step, causing excessive allocation 
          +and runtime. A 1D active-shift specializer briefly improved ho
          +t-loop timings, but it was removed because the core method sho
          +uld remain generic and not be tuned to the target MFIM example
          +.                                                             
      120 +                                                              
      121  ## Validation goals
      122  
      123  For each backend:

● Added a concise “Attempts tried during exploration” section to:

  benchmark/ts_native_trotter_theory.md

  It notes the Taylor-2 diagnostic, naive component-exponential prototype, and
  removed 1D active-shift specializer.

✻ Baked for 38s

※ recap: We’re designing a generic TS-native second-order Trotter method using orbit Liouvillian exponentials over Pauli-orbit graphs. Next, implement or optimize a generic orbit-flow backend beyond the current Krylov baseline.

❯ make a git commit now                                                         
  Committed 3d62611, ran 2 shell commands (ctrl+o to expand)

● Committed the changes.

  3d62611 Add TS-native Trotter prototype

  Included:

  - TrotterTS Krylov orbit-flow prototype
  - generic orbit_edges / orbit_liouvillian
  - evolution tests for TrotterTS
  - slowdown findings markdown
  - TS-native Trotter theory proposal markdown

✻ Worked for 1m 11s

❯ ok lets continue                                                              

● Ready. The code is now cleanly committed with the theory/prototype baseline.
  Next I’d continue from the committed state by investigating a generic cached 
  orbit-graph plan (not a 1D specializer):

  1. build/cache orbit_edges(Ha, q) results per orbit term and TS string,
  2. measure whether repeated Krylov/component calls reuse enough graph data,
  3. only then consider a component-exponential backend again.

  If that sounds right, I’ll start by adding instrumentation/prototype cache
  outside the core API first, so we can verify reduced edge-generation work and
  allocations before changing the implementation.

✻ Worked for 20s

❯ Yes continue                                                                  

● I’ll prototype a generic per-orbit orbit_edges cache outside the package API
  first, then compare it to the committed TrotterTS Krylov path. This will tell
  us whether caching graph generation is actually worth baking into the
  implementation.

● Write(~/.claude/jobs/513a65b3/tmp/cached_orbit_krylov.jl)
  ⎿  Wrote 130 lines to
     ../../../.claude/jobs/513a65b3/tmp/cached_orbit_krylov.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
       8 function Xtot(N)
       9     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
      10 end
     … +142 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I prototyped a generic per-orbit orbit_edges cache outside the package API and
  benchmarked it against committed TrotterTS.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/cached_orbit_krylov.jl

  Results

  ┌────────────────┬───────────────┬──────────────┬────────────┬───────────┐
  │      Case      │   Current     │   Cached     │  Current   │  Cached   │
  │                │   TrotterTS   │  prototype   │   alloc    │   alloc   │
  ├────────────────┼───────────────┼──────────────┼────────────┼───────────┤
  │ N=16, M=10,    │        0.222s │       0.340s │   0.487 GB │  0.570 GB │
  │ steps=10       │               │              │            │           │
  ├────────────────┼───────────────┼──────────────┼────────────┼───────────┤
  │ N=32, M=12,    │        0.189s │       0.327s │   0.430 GB │  0.520 GB │
  │ steps=5        │               │              │            │           │
  └────────────────┴───────────────┴──────────────┴────────────┴───────────┘

  Cache statistics:

  N=16:
    entries = 4656
    hits    = 164852
    misses  = 4656
    edge_time = 0.0017s

  N=32:
    entries = 16213
    hits    = 120025
    misses  = 16213
    edge_time = 0.0095s

  Interpretation

  Caching orbit_edges is not worth adding right now.

  Even though hit rates are high, edge generation itself is not the bottleneck:

  edge_time is tiny relative to total runtime

  The cached version is slower and allocates more because:

  - cache dictionary overhead,
  - storing many small OperatorTS edge objects,
  - extra lookups,
  - still doing expensive Krylov orthogonalization/operator combinations.

  Conclusion

  Do not bake edge caching into the implementation yet.

  The current generic orbit_edges interface is useful for theory/measurement,
  but performance bottleneck is not edge generation. It is likely:

  1. Krylov basis construction/orthogonalization,
  2. repeated trace_product inner products,
  3. repeated Operator arithmetic allocations in _axpy_combination and Arnoldi
  residual updates.

  Next useful direction: instrument inside _orbit_flow_krylov to split
  time/allocations into:

  orbit_liouvillian applications
  trace_product inner products
  operator updates/subtractions
  final axpy combination

  That will show where the Krylov implementation actually loses.

✻ Worked for 1m 26s

❯ yes that is what you should do                                                

● Write(~/.claude/jobs/513a65b3/tmp/krylov_breakdown.jl)
  ⎿  Wrote 133 lines to
     ../../../.claude/jobs/513a65b3/tmp/krylov_breakdown.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +144 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I profiled one full TrotterTS(krylovdim=4) step and split the Krylov
  orbit-flow cost.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/krylov_breakdown.jl

  N=16, M=10, lenO=1024

  One full TS Strang step:

  total = 0.0213s
  alloc = 0.058 GB

  Breakdown:

  ┌──────────────────────────────┬───────────┬───────────────┐
  │             Part             │   Time    │     Calls     │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ orbit_liouvillian            │  0.00809s │            20 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ trace_product inner products │  0.00475s │           100 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ operator updates w - V[i]*h  │  0.00354s │            50 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ truncation                   │  0.00216s │            75 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ final axpy combination       │  0.00139s │ 5 orbit flows │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ norms                        │  0.00091s │            20 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ small dense exp              │ 0.000074s │             5 │
  └──────────────────────────────┴───────────┴───────────────┘

  Public API one-step:

  total = 0.0251s
  alloc = 0.076 GB

  N=32, M=12, lenO=4096

  One full TS Strang step:

  total = 0.127s
  alloc = 0.271 GB

  Breakdown:

  ┌──────────────────────────────┬───────────┬───────────────┐
  │             Part             │   Time    │     Calls     │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ orbit_liouvillian            │   0.0534s │            20 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ trace_product inner products │   0.0291s │           100 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ operator updates w - V[i]*h  │   0.0211s │            50 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ truncation                   │  0.00846s │            75 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ final axpy combination       │  0.00652s │ 5 orbit flows │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ norms                        │  0.00607s │            20 │
  ├──────────────────────────────┼───────────┼───────────────┤
  │ small dense exp              │ 0.000138s │             5 │
  └──────────────────────────────┴───────────┴───────────────┘

  Public API one-step:

  total = 0.126s
  alloc = 0.279 GB

  Interpretation

  The bottleneck is not the small Krylov matrix exponential. It is:

  1. orbit_liouvillian applications
  2. Krylov orthogonalization inner products
  3. allocation-heavy operator linear combinations

  For N=32, approximate time shares:

  orbit_liouvillian: ~42%
  inner products:    ~23%
  operator updates:  ~17%
  truncation:         ~7%
  axpy:               ~5%

  So optimizing only orbit_edges/graph generation is insufficient. Even if orbit
  Liouvillian were free, Krylov orthogonalization and operator arithmetic would
  still be substantial.

  Next obvious engineering target

  The update step:

  w = truncation(w - V[i] * Hk[i,j])

  and final:

  O = O + basis[j] * coeffs[j]

  allocate whole new OperatorTS objects repeatedly.

  A non-mutating implementation is clean but expensive. The next improvement
  should be a small set of internal dictionary-based linear-combination helpers
  for OperatorTS, e.g.:

  linear_combination([(a1, O1), (a2, O2), ...])
  axpy_operator(O, a, V)

  that accumulates directly into one dictionary instead of creating temporary
  operators at each algebraic expression.

  This would target the real hot loop while staying generic.

✻ Worked for 1m 28s

❯ yes                                                                           

● I’ll add generic dictionary-based linear-combination helpers inside
  evolution.jl and use them in the Krylov residual updates and final basis
  combination, then rerun the same breakdown.

● Update(src/evolution.jl)
  ⎿  Added 15 lines, removed 5 lines
      299      return orbit_liouvillian(Ha, O; hbar=hbar)
      300  end
      301  
      302 -function _axpy_combination(basis, coeffs, template)           
      303 -    O = zero(template)                                        
      304 -    for j in eachindex(coeffs)                                
      305 -        iszero(coeffs[j]) || (O = O + basis[j] * coeffs[j])   
      302 +function _operator_linear_combination(template, ops, coeffs; e
          +psilon=1e-16)                                                 
      303 +    d = emptydict(template)                                   
      304 +    for (O, a) in zip(ops, coeffs)                            
      305 +        iszero(a) && continue                                 
      306 +        for (p, c) in zip(O.strings, O.coeffs)                
      307 +            setwith!(+, d, p, a * c)                          
      308 +        end                                                   
      309      end
      307 -    return O                                                  
      310 +    out = typeof(template)(collect(keys(d)), collect(values(d)
          +))                                                            
      311 +    return (eltype(out.coeffs) == ComplexF64) ? cutoff(out, ep
          +silon) : out                                                  
      312  end
      313  
      314 +_operator_axpy(Y, a, X) = _operator_linear_combination(Y, (Y, 
          +X), (one(a), a))                                              
      315 +                                                              
      316 +function _axpy_combination(basis, coeffs, template)           
      317 +    return _operator_linear_combination(template, basis, coeff
          +s)                                                            
      318 +end                                                           
      319 +                                                              
      320  function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
           t, truncation)
      321      krylovdim > 0 || throw(ArgumentError("krylovdim must be po
           sitive"))
      322      β = sqrt(real(trace_product(O', O)))

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      331          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      332          for i in 1:j
      333              Hk[i, j] = trace_product(V[i]', w) / trace_product
           (V[i]', V[i])
      334 -            w = truncation(w - V[i] * Hk[i, j])               
      334 +            w = truncation(_operator_axpy(w, -Hk[i, j], V[i]))
      335          end
      336          nw = sqrt(real(trace_product(w', w)))
      337          if nw < 1e-12

  Ran 2 shell commands (ctrl+o to expand)

● The linear-combination helper did not materially improve the Krylov hot loop.
  It slightly reduced one public-step allocation in one case, but update/axpy
  time did not improve and code clarity got worse. I’m going to revert that
  helper and keep the implementation clean unless a stronger optimization shows
  up.

● Update(src/evolution.jl)
  ⎿  Added 5 lines, removed 15 lines
      299      return orbit_liouvillian(Ha, O; hbar=hbar)
      300  end
      301  
      302 -function _operator_linear_combination(template, ops, coeffs; e
          -psilon=1e-16)                                                 
      303 -    d = emptydict(template)                                   
      304 -    for (O, a) in zip(ops, coeffs)                            
      305 -        iszero(a) && continue                                 
      306 -        for (p, c) in zip(O.strings, O.coeffs)                
      307 -            setwith!(+, d, p, a * c)                          
      308 -        end                                                   
      302 +function _axpy_combination(basis, coeffs, template)           
      303 +    O = zero(template)                                        
      304 +    for j in eachindex(coeffs)                                
      305 +        iszero(coeffs[j]) || (O = O + basis[j] * coeffs[j])   
      306      end
      310 -    out = typeof(template)(collect(keys(d)), collect(values(d)
          -))                                                            
      311 -    return (eltype(out.coeffs) == ComplexF64) ? cutoff(out, ep
          -silon) : out                                                  
      307 +    return O                                                  
      308  end
      309  
      314 -_operator_axpy(Y, a, X) = _operator_linear_combination(Y, (Y, 
          -X), (one(a), a))                                              
      315 -                                                              
      316 -function _axpy_combination(basis, coeffs, template)           
      317 -    return _operator_linear_combination(template, basis, coeff
          -s)                                                            
      318 -end                                                           
      319 -                                                              
      310  function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
           t, truncation)
      311      krylovdim > 0 || throw(ArgumentError("krylovdim must be po
           sitive"))
      312      β = sqrt(real(trace_product(O', O)))

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      321          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      322          for i in 1:j
      323              Hk[i, j] = trace_product(V[i]', w) / trace_product
           (V[i]', V[i])
      324 -            w = truncation(_operator_axpy(w, -Hk[i, j], V[i]))
      324 +            w = truncation(w - V[i] * Hk[i, j])               
      325          end
      326          nw = sqrt(real(trace_product(w', w)))
      327          if nw < 1e-12
● I tried the generic dictionary-based linear-combination optimization and
  measured it.

  What I changed experimentally

  Replaced Krylov residual updates:

  w = truncation(w - V[i] * Hk[i,j])

  with a dictionary-accumulating helper to avoid constructing intermediate
  scaled/subtracted operators.

  Also used the helper for final Krylov basis combination.

  Breakdown after change

  N=16, M=10

  Before:

  counted step total: 0.0213s
  alloc: 0.058 GB
  update: 0.00354s
  axpy:   0.00139s
  public step alloc: 0.076 GB

  After helper:

  counted step total: 0.0218s
  alloc: 0.058 GB
  update: 0.00367s
  axpy:   0.00147s
  public step alloc: 0.080 GB

  Worse/slightly neutral.

  N=32, M=12

  Before:

  counted step total: 0.1270s
  alloc: 0.271 GB
  update: 0.0211s
  axpy:   0.00652s
  public step alloc: 0.279 GB

  After helper:

  counted step total: 0.1252s
  alloc: 0.271 GB
  update: 0.0232s
  axpy:   0.00810s
  public step alloc: 0.266 GB

  Slight total noise-level improvement and small public allocation reduction,
  but the targeted update/axpy operations got slower.

  Decision

  I reverted the helper.

  Reason: it does not clearly improve the hot loop and makes the code less
  clear. As you noted, preallocation/low-level engineering is only worth it when
  it is truly buying something in a hot loop. This one didn’t.

  Current bottleneck remains

  From the clean breakdown:

  N=32 one TrotterTS step:
    orbit_liouvillian: ~42%
    inner products:    ~23%
    operator updates:  ~17%
    truncation:         ~7%

  The next meaningful direction is probably not generic linear-combination
  syntax, but either:

  1. reduce Krylov orthogonalization cost, or
  2. find a cheaper orbit-exponential backend than global Krylov, with component
  plans/caching done carefully.

✻ Cooked for 2m 32s

❯ Attempt to reduce the allocations once more                                   

● I’ll make one allocation-focused change that preserves the generic graph
  interface: keep orbit_edges for inspection/proposals, but stop using it inside
  orbit_liouvillian because that creates many small temporary OperatorTS edge
  objects. The Liouvillian hot path should use the existing generic TS
  commutator directly.

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 10 lines
      278  
      279  function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
           gth::Int=1000)
      280      checklength(A, B)
      281 -    d = emptydict(A)                                          
      282 -    for (q, cq) in zip(B.strings, B.coeffs)                   
      283 -        edges = orbit_edges(A, q; epsilon=epsilon, maxlength=m
          -axlength)                                                     
      284 -        for (p, c) in zip(edges.strings, edges.coeffs)        
      285 -            setwith!(+, d, p, cq * c)                         
      286 -        end                                                   
      287 -    end                                                       
      288 -    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      289 -    o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) 
          -: o                                                           
      290 -    return 1im * o / hbar                                     
      281 +    return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
          +axlength) / hbar                                              
      282  end
      283  
      284  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)

  Ran 2 shell commands (ctrl+o to expand)

● Good allocation reduction this time.

  Change made

  Kept the generic orbit_edges API for graph inspection, but stopped using it
  inside the hot orbit_liouvillian path.

  Before:

  orbit_liouvillian(A, B)
      for q in B
          edges = orbit_edges(A, q)  # many small temporary OperatorTS objects
          accumulate edges

  After:

  orbit_liouvillian(A, B) =
      1im * commutator(A, B) / hbar

  This preserves the public graph interface while avoiding per-string temporary
  edge operators in Krylov.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

  Krylov breakdown improvement

  N=16, M=10

  Before:

  counted step total: 0.0213s
  alloc: 0.058 GB
  liouvillian: 0.00809s
  public step: 0.0251s, 0.076 GB

  After:

  counted step total: 0.0183s
  alloc: 0.041 GB
  liouvillian: 0.00471s
  public step: 0.0214s, 0.060 GB

  N=32, M=12

  Before:

  counted step total: 0.1270s
  alloc: 0.271 GB
  liouvillian: 0.0534s
  public step: 0.1257s, 0.279 GB

  After:

  counted step total: 0.0954s
  alloc: 0.197 GB
  liouvillian: 0.0283s
  public step: 0.1013s, 0.205 GB

  Interpretation

  This was worth keeping:

  - N=32 one-step allocation dropped by ~27%.
  - N=32 one-step time dropped by ~20–25%.
  - The bottleneck shifted: inner products now dominate slightly more than
  orbit_liouvillian.

  Current N=32 counted breakdown:

  orbit_liouvillian: 0.0283s
  inner products:    0.0308s
  operator updates:  0.0166s
  truncation:        0.0057s
  axpy:              0.0056s

  So next allocation/time target is the trace_product inner products and
  orthogonalization structure, not graph edge generation.

✻ Worked for 1m 21s

❯ go ahead                                                                      

  Read 2 files (ctrl+o to expand)

● I’ll target the inner-product allocations now. The current Krylov code calls
  trace_product(V[i]', w), which allocates an adjoint operator every time. I’ll
  add an internal adjoint-aware trace helper and use it in Krylov, then measure.

● Update(src/evolution.jl)
  ⎿  Added 30 lines
      299      return orbit_liouvillian(Ha, O; hbar=hbar)
      300  end
      301  
      302 +function _trace_product_adjoint_left(o1::Operator{<:PauliStrin
          +gTS}, o2::Operator{<:PauliStringTS})                          
      303 +    checklength(o1, o2)                                       
      304 +    Ls = qubitsize(o1)                                        
      305 +    Ps = periodicflags(o1)                                    
      306 +    tr = zero(promote_type(scalartype(o1), scalartype(o2)))   
      307 +                                                              
      308 +    d = emptydict(o2)                                         
      309 +    for (p2, c2) in zip(o2.strings, o2.coeffs)                
      310 +        insert!(d, p2, c2)                                    
      311 +    end                                                       
      312 +                                                              
      313 +    for (p1, c1) in zip(o1.strings, o1.coeffs)                
      314 +        c2 = get(d, p1, nothing)                              
      315 +        isnothing(c2) && continue                             
      316 +        rep1 = representative(p1)                             
      317 +        p, k = prod(rep1, rep1)                               
      318 +        s = 1 - ((ycount(p1) & 1) << 1)                       
      319 +        f = s * conj(c1) * c2 * k                             
      320 +        for shift_tuple in all_shifts(Ls, Ps)                 
      321 +            shifted = shift(rep1, Ls, Ps, shift_tuple)        
      322 +            if shifted == rep1                                
      323 +                tr += f                                       
      324 +            end                                               
      325 +        end                                                   
      326 +    end                                                       
      327 +    scale = 2.0^Base.prod(Ls)                                 
      328 +    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) i
          +f p)                                                          
      329 +    return tr * scale * num_translations                      
      330 +end                                                           
      331 +                                                              
      332  function _axpy_combination(basis, coeffs, template)
      333      O = zero(template)
      334      for j in eachindex(coeffs)

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      339  
      340  function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
           t, truncation)
      341      krylovdim > 0 || throw(ArgumentError("krylovdim must be po
           sitive"))
      342 -    β = sqrt(real(trace_product(O', O)))                      
      342 +    β = sqrt(real(_trace_product_adjoint_left(O, O)))         
      343      β == 0 && return O
      344  
      345      V = typeof(O)[]

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      350      for j in 1:krylovdim
      351          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      352          for i in 1:j
      353 -            Hk[i, j] = trace_product(V[i]', w) / trace_product
          -(V[i]', V[i])                                                 
      353 +            Hk[i, j] = _trace_product_adjoint_left(V[i], w) / 
          +_trace_product_adjoint_left(V[i], V[i])                       
      354              w = truncation(w - V[i] * Hk[i, j])
      355          end
      356          nw = sqrt(real(trace_product(w', w)))

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      353              Hk[i, j] = _trace_product_adjoint_left(V[i], w) / 
           _trace_product_adjoint_left(V[i], V[i])
      354              w = truncation(w - V[i] * Hk[i, j])
      355          end
      356 -        nw = sqrt(real(trace_product(w', w)))                 
      356 +        nw = sqrt(real(_trace_product_adjoint_left(w, w)))    
      357          if nw < 1e-12
      358              actual_dim = j
      359              break

  Ran 2 shell commands (ctrl+o to expand)

● I attempted another allocation reduction by removing adjoint allocations in
  Krylov inner products.

  Change tried

  Krylov originally used:

  trace_product(V[i]', w)
  trace_product(V[i]', V[i])
  trace_product(w', w)

  Each ' allocates an adjoint operator.

  I added an internal helper:

  _trace_product_adjoint_left(o1, o2)

  that computes the equivalent of:

  trace_product(o1', o2)

  without constructing o1'.

  Then I replaced the Krylov uses.

  Result

  Tests still pass:

  evolve | 8 passed

  But performance did not improve.

  Before this change, after the previous Liouvillian allocation fix

  N=16:
    counted step total = 0.0183s
    alloc = 0.041 GB
    inner = 0.00496s
    public step = 0.0214s, 0.060 GB

  N=32:
    counted step total = 0.0954s
    alloc = 0.197 GB
    inner = 0.0308s
    public step = 0.1013s, 0.205 GB

  After adjoint-free trace helper

  N=16:
    counted step total = 0.0193s
    alloc = 0.041 GB
    inner = 0.00518s
    public step = 0.0210s, 0.065 GB

  N=32:
    counted step total = 0.0989s
    alloc = 0.197 GB
    inner = 0.0304s
    public step = 0.1027s, 0.199 GB

  Allocations are essentially unchanged. Timing is noise/slightly worse.

  Interpretation

  The adjoint allocation was not a meaningful part of the total allocation. The
  helper duplicates trace_product logic and does not buy enough to justify the
  code complexity.

  I should revert this helper unless you want to keep it for semantic clarity.
  My recommendation is to remove it and return to:

  trace_product(V[i]', w)

  The only allocation reduction that clearly helped was the previous one:

  orbit_liouvillian uses commutator directly instead of orbit_edges per q

  That one should stay.

✻ Cogitated for 2m 31s

❯ Revert. Actually profile and see what the bottlenecks are: are they mainly in 
  the main orbit splitting algorithm, or in the krylov method approximation,    
  or due to heavy allocation in hot loops.                                      

● Update(src/evolution.jl)
  ⎿  Removed 30 lines
      299      return orbit_liouvillian(Ha, O; hbar=hbar)
      300  end
      301  
      302 -function _trace_product_adjoint_left(o1::Operator{<:PauliStrin
          -gTS}, o2::Operator{<:PauliStringTS})                          
      303 -    checklength(o1, o2)                                       
      304 -    Ls = qubitsize(o1)                                        
      305 -    Ps = periodicflags(o1)                                    
      306 -    tr = zero(promote_type(scalartype(o1), scalartype(o2)))   
      307 -                                                              
      308 -    d = emptydict(o2)                                         
      309 -    for (p2, c2) in zip(o2.strings, o2.coeffs)                
      310 -        insert!(d, p2, c2)                                    
      311 -    end                                                       
      312 -                                                              
      313 -    for (p1, c1) in zip(o1.strings, o1.coeffs)                
      314 -        c2 = get(d, p1, nothing)                              
      315 -        isnothing(c2) && continue                             
      316 -        rep1 = representative(p1)                             
      317 -        p, k = prod(rep1, rep1)                               
      318 -        s = 1 - ((ycount(p1) & 1) << 1)                       
      319 -        f = s * conj(c1) * c2 * k                             
      320 -        for shift_tuple in all_shifts(Ls, Ps)                 
      321 -            shifted = shift(rep1, Ls, Ps, shift_tuple)        
      322 -            if shifted == rep1                                
      323 -                tr += f                                       
      324 -            end                                               
      325 -        end                                                   
      326 -    end                                                       
      327 -    scale = 2.0^Base.prod(Ls)                                 
      328 -    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) i
          -f p)                                                          
      329 -    return tr * scale * num_translations                      
      330 -end                                                           
      331 -                                                              
      302  function _axpy_combination(basis, coeffs, template)
      303      O = zero(template)
      304      for j in eachindex(coeffs)

● Agreed. I’ll revert the adjoint-free trace helper, then profile with Julia’s
  allocation/time tooling at a coarser level to attribute cost between: orbit
  splitting orchestration, orbit Liouvillian applications, Krylov
  orthogonalization, and allocation-heavy operator algebra.

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      309  
      310  function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
           Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
           t, truncation)
      311      krylovdim > 0 || throw(ArgumentError("krylovdim must be po
           sitive"))
      312 -    β = sqrt(real(_trace_product_adjoint_left(O, O)))         
      312 +    β = sqrt(real(trace_product(O', O)))                      
      313      β == 0 && return O
      314  
      315      V = typeof(O)[]

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      320      for j in 1:krylovdim
      321          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      322          for i in 1:j
      323 -            Hk[i, j] = _trace_product_adjoint_left(V[i], w) / 
          -_trace_product_adjoint_left(V[i], V[i])                       
      323 +            Hk[i, j] = trace_product(V[i]', w) / trace_product
          +(V[i]', V[i])                                                 
      324              w = truncation(w - V[i] * Hk[i, j])
      325          end
      326          nw = sqrt(real(_trace_product_adjoint_left(w, w)))

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      323              Hk[i, j] = trace_product(V[i]', w) / trace_product
           (V[i]', V[i])
      324              w = truncation(w - V[i] * Hk[i, j])
      325          end
      326 -        nw = sqrt(real(_trace_product_adjoint_left(w, w)))    
      326 +        nw = sqrt(real(trace_product(w', w)))                 
      327          if nw < 1e-12
      328              actual_dim = j
      329              break

● Write(~/.claude/jobs/513a65b3/tmp/profile_krylov_allocs.jl)
  ⎿  Wrote 31 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_krylov_allocs.jl
       1 using PauliStrings
       2 using Profile
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +25 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)
● Reverted the adjoint-free trace helper and profiled the current Krylov
  TrotterTS path.

  Current important code state

  orbit_liouvillian now uses the generic TS commutator directly:

  orbit_liouvillian(A, B; hbar=1, ...) =
      1im * commutator(A, B; ...) / hbar

  This avoids the earlier allocation-heavy implementation that called
  orbit_edges separately for every input string.

  orbit_edges remains available for graph inspection/theory, but is no longer
  used in the hot Krylov Liouvillian path.

  One-step Krylov breakdown

  Measured one full TrotterTS(krylovdim=4) Strang step on an evolved state.

  N=16, M=10, lenO=1024

  counted_step total = 0.0185s
  alloc              = 0.041 GB
  public_step total  = 0.0211s
  public_step alloc  = 0.060 GB

  Breakdown:

  ┌───────────────────────────┬───────────┬───────┐
  │         Component         │   Time    │ Calls │
  ├───────────────────────────┼───────────┼───────┤
  │ orbit_liouvillian         │  0.00474s │    20 │
  ├───────────────────────────┼───────────┼───────┤
  │ Krylov inner products     │  0.00495s │   100 │
  ├───────────────────────────┼───────────┼───────┤
  │ operator residual updates │  0.00374s │    50 │
  ├───────────────────────────┼───────────┼───────┤
  │ truncation                │  0.00218s │    75 │
  ├───────────────────────────┼───────────┼───────┤
  │ final axpy combination    │  0.00143s │     5 │
  ├───────────────────────────┼───────────┼───────┤
  │ norms                     │  0.00099s │    20 │
  ├───────────────────────────┼───────────┼───────┤
  │ small dense exp           │ 0.000044s │     5 │
  └───────────────────────────┴───────────┴───────┘

  N=32, M=12, lenO=4096

  counted_step total = 0.103s
  alloc              = 0.197 GB
  public_step total  = 0.100s
  public_step alloc  = 0.205 GB

  Breakdown:

  ┌───────────────────────────┬──────────┬───────┐
  │         Component         │   Time   │ Calls │
  ├───────────────────────────┼──────────┼───────┤
  │ Krylov inner products     │  0.0334s │   100 │
  ├───────────────────────────┼──────────┼───────┤
  │ orbit_liouvillian         │  0.0295s │    20 │
  ├───────────────────────────┼──────────┼───────┤
  │ operator residual updates │  0.0178s │    50 │
  ├───────────────────────────┼──────────┼───────┤
  │ truncation                │ 0.00771s │    75 │
  ├───────────────────────────┼──────────┼───────┤
  │ norms                     │ 0.00624s │    20 │
  ├───────────────────────────┼──────────┼───────┤
  │ final axpy combination    │ 0.00622s │     5 │
  ├───────────────────────────┼──────────┼───────┤
  │ small dense exp           │ 0.00010s │     5 │
  └───────────────────────────┴──────────┴───────┘

  Sampling profile confirms the same picture

  Hot functions in the profile are mainly:

  - orbit_liouvillian
  - commutator
  - TS binary_kernel
  - trace_product
  - trim / partialsortperm
  - operator scalar multiplication/subtraction in Krylov residual updates

  The small Krylov dense matrix exponential is negligible.

  Answer to the bottleneck question

  Not mainly the orbit splitting algorithm

  The orbit-level Strang orchestration itself is cheap:

  5 orbit factors per step for m=3

  That part is not the issue.

  Mainly the Krylov approximation machinery

  With krylovdim=4, one Strang step performs roughly:

  5 orbit flows
  × 4 Krylov Liouvillian applications
  = 20 orbit_liouvillian calls

  plus:
  100 trace_product inner products
  50 residual operator updates
  75 truncations

  So the current cost is mostly Krylov machinery, not the abstract orbit
  splitting.

  Allocations are significant and in hot loops

  At N=32, one step allocates about:

  ~0.2 GB

  The allocation-heavy parts are:

  - repeated OperatorTS temporaries in Krylov residual updates,
  - trace_product(V[i]', w) and related inner products,
  - orbit_liouvillian returning fresh operators,
  - truncations after each Krylov operation.

  The attempted adjoint-free trace helper did not reduce allocations
  meaningfully, so I reverted it.

  Most useful conclusion

  The current TS-native method is mathematically clean but the naive Krylov
  backend is expensive.

  The bottleneck is not:

  orbit splitting

  It is:

  matrix-free Krylov over sparse Pauli-orbit operators, implemented with
  high-level OperatorTS arithmetic

  The next serious improvement should target the Krylov backend as a whole,
  likely by making Krylov operate on a coefficient vector over an
  explicit/cached orbit graph basis, rather than repeatedly constructing full
  OperatorTS intermediates for each inner product/update.

✻ Sautéed for 2m 55s

❯ try that                                                                      

● I’ll prototype a vector-space Krylov backend outside the package first: sparse
  coefficient dictionaries over TS Pauli-orbit keys, orbit_edges as the matrix
  action, and weighted Hilbert-Schmidt inner products. Then I’ll compare
  one-step and short-run timing/allocations against the current OperatorTS
  Krylov path.

● Write(~/.claude/jobs/513a65b3/tmp/vector_krylov_prototype.jl)
  ⎿  Wrote 173 lines to
     ../../../.claude/jobs/513a65b3/tmp/vector_krylov_prototype.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
       8 function Xtot(N)
       9     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
      10 end
     … +185 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I prototyped the vector-space Krylov backend outside the package.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/vector_krylov_prototype.jl

  Prototype idea

  Instead of representing Krylov vectors as OperatorTS, represent them as:

  Dict{PauliStringTS, ComplexF64}

  and implement:

  - liouv_vec using orbit_edges,
  - weighted inner products directly over coefficient dictionaries,
  - axpy! directly in coefficient dictionaries,
  - convert to OperatorTS only for trim.

  This was meant to reduce full OperatorTS temporaries in Krylov.

  Results

  N=16, M=10

  ┌──────────────────────────┬─────────┬──────────┐
  │          Method          │  Time   │  Alloc   │
  ├──────────────────────────┼─────────┼──────────┤
  │ Vector Krylov prototype  │ 0.0445s │ 0.049 GB │
  ├──────────────────────────┼─────────┼──────────┤
  │ Current public TrotterTS │ 0.0180s │ 0.169 GB │
  └──────────────────────────┴─────────┴──────────┘

  Vector version reduced allocation but was ~2.5x slower.

  Stats:

  edge_time = 0.0021s
  edge_calls = 4075
  dot_time = 0.0056s
  dot_calls = 125
  trim_time = 0.0038s
  cache_entries = 4075
  metric_entries = 1927

  N=32, M=12

  ┌──────────────────────────┬────────┬──────────┐
  │          Method          │  Time  │  Alloc   │
  ├──────────────────────────┼────────┼──────────┤
  │ Vector Krylov prototype  │ 0.194s │ 0.221 GB │
  ├──────────────────────────┼────────┼──────────┤
  │ Current public TrotterTS │ 0.100s │ 0.232 GB │
  └──────────────────────────┴────────┴──────────┘

  Allocation is similar/slightly better, but runtime is ~2x worse.

  Stats:

  edge_time = 0.0084s
  edge_calls = 16111
  dot_time = 0.0193s
  dot_calls = 125
  trim_time = 0.0137s
  cache_entries = 16111
  metric_entries = 7451

  Interpretation

  The vector-space direction reduces some allocations, but the naive dictionary
  implementation is too slow.

  The cost moves into:

  - Julia Dict lookups/updates,
  - repeated conversion to OperatorTS for trimming,
  - orbit_edges calls for every key,
  - metric-cache lookup per inner-product entry.

  So a vector-space backend might only be viable if implemented with compact
  indexed arrays over a fixed/cached basis, not arbitrary dictionaries.

  Conclusion

  Do not implement this prototype in the package.

  Current OperatorTS Krylov is faster despite higher allocations because its
  existing dictionary kernels are better optimized than this ad hoc
  coefficient-dict prototype.

  A real vector backend would need:

  basis indexing: PauliStringTS -> Int
  coeff vectors: Vector{ComplexF64}
  sparse matrix/edge lists: Vector{(target_index, weight)}

  i.e. a proper cached orbit graph plan, not plain dictionaries.

✻ Crunched for 1m 56s

❯ I don't think we should be reimplementing primitives already in the package.  
  We can assume that these are already best effort, and those optimizations lie 
  outside the scope of this PR. Let us revisit what we are doing                

● Agreed. We’ve drifted into reimplementing lower-level primitives
  (trace_product, linear combinations, coefficient-vector kernels) that the
  package already owns. That is outside scope and risks making the PR both less
  clear and less maintainable.

  Let’s reset the scope.

  What this PR should be about

  The original task is:

  ▎ Design a second-order Trotter method that exploits translation invariance 
  ▎ and can outperform RK4 in the TS case.

  The core contribution should be algorithmic:

  1. diagnose why current TS Trotter is slow,
  2. introduce the correct TS-native orbit-level formulation,
  3. add clean primitives/prototype implementation that express that
  formulation,
  4. document findings and theory,
  5. avoid target-specific or low-level overengineering.

  What we have established

  1. Current TS Trotter slowdown

  The slowdown is not mainly resum(H) allocation. It is that current TS Trotter
  turns m TS Hamiltonian orbits into O(Nm) translated ordinary gates.

  For MFIM:

  m = 3
  expanded Strang gates = 6N - 1

  At N=64:

  383 ordinary gate updates per step

  This is the wrong granularity.

  2. Correct TS-level structure

  The natural TS-native split is:

  H = \sum_a H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_r T_r P_a T_r^\dagger.

  Then use Strang over orbit Liouvillians:

  e^{dt/2 L_1}
  e^{dt/2 L_2}
  ...
  e^{dt L_m}
  ...
  e^{dt/2 L_1}

  where:

  L_a(O) = i[H_a^{orb}, O] / \hbar.

  This reduces product-formula factors from:

  O(Nm)

  to:

  O(m)

  3. Proper orbit-flow primitive

  The core orbit flow is:

  O \mapsto e^{τ L_a} O.

  This is not the same as RK4 if treated as an exponential action. It is a real
  product formula over grouped TS orbit Hamiltonians.

  4. Taylor-2 was diagnostic only

  Taylor-2 was fast-ish but conceptually too close to RK/RK2-style commutator
  integration. It should not be in the final code.

  5. Krylov is conceptually right but slow

  TrotterTS(krylovdim=4) is clean and generic, but currently slower than RK4
  because it invokes many OperatorTS operations inside Arnoldi.

  That is acceptable for a prototype if we frame it honestly:

  the PR establishes the TS-native orbit Trotter formulation and a
  correctness-oriented Krylov implementation;
  performance optimization of the underlying orbit exponential backend is future
  work.

  6. Graph primitives are useful

  orbit_edges and orbit_liouvillian are useful conceptual primitives:

  orbit_edges(Ha, q)
  orbit_liouvillian(Ha, O)

  They expose the Pauli-orbit graph and allow validation/measurement.

  But the hot implementation should not reimplement existing low-level package
  operations. For now, orbit_liouvillian can simply wrap the existing optimized
  commutator:

  1im * commutator(A, B) / hbar

  and orbit_edges can remain an inspection/proposal primitive.

  What should stay in the PR

  I think these should stay:

  Code

  - TrotterTS(; order=2, krylovdim=4)
  - _orbit_terms
  - _trotterts_step
  - _orbit_flow_krylov
  - orbit_liouvillian as a generic wrapper around commutator
  - maybe orbit_edges if we want to expose the graph rule

  Tests

  - small exact comparison for TrotterTS,
  - output remains OperatorTS.

  Docs/notes

  - benchmark/ts_trotter_slowdown_findings.md
  - benchmark/ts_native_trotter_theory.md

  What should not stay

  Remove or avoid:

  - Taylor-2 implementation,
  - 1D active-shift specializer,
  - local component exponential prototype,
  - coefficient-dict vector Krylov prototype,
  - custom trace-product helpers,
  - custom linear-combination helpers.

  We already removed most of these. Good.

  Revised positioning

  The PR is not “final optimized TS Trotter beats RK4 everywhere.”

  The PR is:

  ▎ Introduce a TS-native orbit-level Trotter formulation and prototype 
  ▎ implementation. The previous TS Trotter used the wrong granularity (resum(H)
  ▎ → translated ordinary gates). This PR establishes the correct orbit 
  ▎ Liouvillian split and documents why further performance work should target 
  ▎ the orbit exponential backend.

  That is much more honest and maintainable.

  Next cleanup step

  I’d now do a final code review for scope:

  1. Confirm only clean generic code remains.
  2. Run tests.
  3. Maybe adjust docs/comments to say TrotterTS is experimental/prototype if
  needed.
  4. Commit the cleanup separately.

  Then later, if desired, we can open a separate performance PR for a true
  optimized orbit exponential backend.

✻ Worked for 41s

❯ No this is not what I meant. The end goal is still to complete the TrotterTS, 
  possibly using orbit splitting. I meant, that if slowdowns are occuring due   
  to multiple intermediate OperatorTS creations, then we should focus on        
  reducing these                                                                

● Understood. We should keep pushing toward a complete TrotterTS; the question
  is whether the current Krylov backend is slow because it creates too many
  intermediate OperatorTS objects. I’ll target that specifically while still
  using existing package primitives.

  First low-risk fix: cache Krylov basis adjoints and norms so we don’t
  repeatedly allocate V[i]' and recompute trace_product(V[i]', V[i]) in the
  inner loop.

  Read 1 file (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 11 lines, removed 3 lines
      313      β == 0 && return O
      314  
      315      V = typeof(O)[]
      316 -    push!(V, O / β)                                           
      316 +    Vdag = typeof(O)[]                                        
      317 +    Vnorms = ComplexF64[]                                     
      318 +    v1 = O / β                                                
      319 +    push!(V, v1)                                              
      320 +    push!(Vdag, v1')                                          
      321 +    push!(Vnorms, trace_product(Vdag[1], V[1]))               
      322      Hk = zeros(ComplexF64, krylovdim + 1, krylovdim)
      323      actual_dim = krylovdim
      324  
      325      for j in 1:krylovdim
      326          w = truncation(_orbit_liouvillian(Ha, V[j], hbar))
      327          for i in 1:j
      323 -            Hk[i, j] = trace_product(V[i]', w) / trace_product
          -(V[i]', V[i])                                                 
      328 +            Hk[i, j] = trace_product(Vdag[i], w) / Vnorms[i]  
      329              w = truncation(w - V[i] * Hk[i, j])
      330          end
      331          nw = sqrt(real(trace_product(w', w)))
     ...
      334              break
      335          end
      336          Hk[j + 1, j] = nw
      332 -        push!(V, w / nw)                                      
      337 +        v = w / nw                                            
      338 +        push!(V, v)                                           
      339 +        push!(Vdag, v')                                       
      340 +        push!(Vnorms, trace_product(Vdag[end], V[end]))       
      341      end
      342  
      343      F = exp(dt * Hk[1:actual_dim, 1:actual_dim])

  Ran 1 shell command (ctrl+o to expand)

● Tried a targeted reduction of intermediate OperatorTS creations in Krylov
  inner products by caching Krylov basis adjoints and norms.

  Change

  Inside _orbit_flow_krylov, instead of repeatedly doing:

  trace_product(V[i]', w)
  trace_product(V[i]', V[i])

  for every inner product, it now stores:

  Vdag[i] = V[i]'
  Vnorms[i] = trace_product(Vdag[i], V[i])

  when each Krylov vector is created.

  Then the inner loop uses:

  trace_product(Vdag[i], w) / Vnorms[i]

  This avoids repeated adjoint construction and denominator recomputation.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

  Profiling result

  N=16, M=10

  Before:

  counted total: 0.01847s
  alloc:         0.041 GB
  public step:   0.02111s, 0.060 GB
  inner:         0.00495s

  After:

  counted total: 0.01877s
  alloc:         0.041 GB
  public step:   0.02232s, 0.056 GB
  inner:         0.00507s

  N=32, M=12

  Before:

  counted total: 0.10296s
  alloc:         0.197 GB
  public step:   0.10006s, 0.205 GB
  inner:         0.03336s

  After:

  counted total: 0.10234s
  alloc:         0.197 GB
  public step:   0.09617s, 0.185 GB
  inner:         0.03242s

  Interpretation

  This one is worth keeping, but only mildly:

  - It reduces public-step allocations for N=32 by ~10%.
  - It gives a small public-step time improvement in that case.
  - It does not change the counted inner-loop allocation much, so it is not a
  breakthrough.

  Still, unlike the previous helper, it is simple and directly reduces repeated
  intermediate adjoint creation in a hot loop. It’s a reasonable cleanup.

  Current bottleneck after this

  At N=32, one Krylov TS-Trotter step is still roughly:

  liouvillian: 0.0315s
  inner:       0.0324s
  updates:     0.0159s
  truncation:  0.0062s
  axpy:        0.0069s
  norms:       0.0072s

  So the bottleneck is split between:

  1. orbit Liouvillian applications,
  2. Krylov orthogonalization inner products,
  3. residual updates.

  The orbit splitting itself remains cheap; the Krylov backend is the expensive
  part.

✻ Worked for 1m 39s

❯ Ok. How are we applying the orbit Louvillian                                  

● Right now the actual hot path applies the orbit Liouvillian through the
  existing TS commutator machinery.

  The call chain is:

  TrotterTS evolve
    -> _trotterts_step(...)
    -> _orbit_flow(...)
    -> _orbit_flow_krylov(...)
    -> _orbit_liouvillian(Ha, V[j], hbar)
    -> orbit_liouvillian(Ha, V[j]; hbar=hbar)
    -> 1im * commutator(Ha, V[j]) / hbar

  Concretely, in src/evolution.jl:

  function _orbit_liouvillian(Ha::Operator{<:PauliStringTS},
                              O::Operator{<:PauliStringTS},
                              hbar::Real)
      return orbit_liouvillian(Ha, O; hbar=hbar)
  end

  and in src/translation_symmetry.jl:

  function orbit_liouvillian(A::Operator{<:PauliStringTS},
                             B::Operator{<:PauliStringTS};
                             hbar::Real=1,
                             epsilon::Real=0,
                             maxlength::Int=1000)
      checklength(A, B)
      return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
  end

  So for the current implementation:

  L_A(B) = \frac{i}{\hbar}[A, B]

  is just delegated to the package’s existing commutator(A, B).

  What commutator(A, B) does for TS operators

  For:

  A::Operator{<:PauliStringTS}
  B::Operator{<:PauliStringTS}

  the method dispatches to the TS binary_kernel in src/translation_symmetry.jl.

  That kernel does roughly:

  for (p1, c1) in A
      rep1 = representative(p1)

      for (p2, c2) in B
          rep2 = representative(p2)

          for s in all_shifts(paulistringtype(A))
              shifted = shift(rep2, Ls, Ps, s)
              p, k = commutator(rep1, shifted)

              if k != 0
                  setwith!(+, d, PauliStringTS{Ls,Ps}(p), c1 * c2 * k)
              end
          end
      end
  end

  Then it constructs:

  OperatorTS{Ls,Ps}(...)

  from the accumulated dictionary and applies cutoff.

  So mathematically it implements:

  [A, B]
  =
  \sum_{p_1 \in A}
  \sum_{p_2 \in B}
  \sum_{\delta \in G}
  c_1 c_2 [p_1, T_\delta p_2 T_\delta^\dagger]

  and folds each result back to a canonical TS representative.

  Important detail

  For TrotterTS, each orbit term is a single representative:

  Ha = typeof(H)([H.strings[j]], [H.coeffs[j]])

  so the hot Liouvillian call is usually:

  commutator(Ha, O)

  where:

  length(Ha) = 1
  length(O) = current number of TS strings

  Thus each Liouvillian application loops over:

  length(O) × number_of_translations

  relative shifts, then folds outputs into TS strings.

  What orbit_edges does

  We also have:

  orbit_edges(Ha, q)

  which exposes the same graph rule for a single input orbit q.

  But it is not used in the hot TrotterTS path anymore, because using it inside
  orbit_liouvillian created many tiny temporary OperatorTS objects and was
  slower.

  So currently:

  - orbit_edges: graph inspection / future component planning
  - orbit_liouvillian: hot path, delegates to commutator

  Summary

  The orbit Liouvillian is currently applied as:

  orbit_liouvillian(Ha, O) = 1im * commutator(Ha, O) / hbar

  and the TS commutator computes the translation-orbit action by scanning all
  relative shifts and accumulating folded TS representatives.

✻ Brewed for 27s

❯ Ok. If the commutator is decomposing the OperatorTS into all the pauli        
  strings, isn't this basically the same issue as resum? If we are to get a TS  
  native update, then the commutator must also remain TS native? Please correct 
  me if I am wrong                                                              

● You’re not wrong to worry, but there is an important distinction:

  resum(H) expands into all translated ordinary Pauli strings as state/gates.
  TS commutator loops over relative translations but accumulates immediately
  back into OperatorTS.

  So the current TS commutator is TS-native in representation, but not fully
  TS-native in algorithmic cost.

  What the current TS commutator does

  For:

  commutator(A::OperatorTS, B::OperatorTS)

  it does not literally build:

  resum(A)
  resum(B)

  and then compute the full ordinary commutator.

  Instead, it computes in the quotient/orbit basis:

  for representative p1 in A
      for representative p2 in B
          for relative shift δ
              compute [p1, shift(p2, δ)]
              fold result back to PauliStringTS

  So the output remains compressed as OperatorTS throughout.

  That is why RK4 in the TS case is much better than doing RK4 on resum(H), 
  resum(O).

  But yes: it still scans all translations

  The cost still includes:

  #representatives(A) × #representatives(B) × #translations

  For Ha as one orbit term:

  length(Ha) = 1

  so one orbit Liouvillian call costs roughly:

  length(O) × N

  relative-shift checks.

  This is not as bad as full resum expansion of both sides, which would be more
  like:

  (length(A)×N) × (length(B)×N)

  with ordinary strings, but it still pays an N factor.

  So the precise diagnosis is

  Current commutator(OperatorTS, OperatorTS) is:

  representation-native: yes
  cost-native: only partially

  It avoids materializing the full expanded operator, but it still enumerates
  every relative translation.

  A more fully TS-native Pauli-orbit update would avoid scanning translations
  that cannot contribute, and ideally work directly with the sparse orbit graph:

  [Q] -> small set of [R] edges

  rather than:

  for every δ in translations
      check whether [P_δ, Q] is nonzero

  Why this matters for TrotterTS

  Our current TrotterTS applies orbit flows with Krylov. Krylov repeatedly
  calls:

  orbit_liouvillian(Ha, V[j])

  and that currently calls:

  commutator(Ha, V[j])

  So even though the Trotter splitting is orbit-level, each Liouvillian
  application still does a full relative-shift scan.

  That is the remaining gap.

  What “fully TS-native” should mean here

  For one orbit generator [P] and one input orbit [Q], the update rule should
  be:

  L_P([Q]) = \sum_{R} a_{R,Q} [R]

  where edges Q -> R are generated only from active relative placements where P
  overlaps/anticommutes with Q.

  So the graph degree should scale with relative support, not total system size:

  current: O(N) relative shifts per Q
  desired: O(number of active relative overlaps)

  Our measurements showed that for the tested states at N=64, only around:

  ~6–7% of shifts actually produced nonzero commutators

  and folded graph degree was around:

  ~4 outgoing edges per input orbit per Hamiltonian orbit

  So a true TS-native commutator/orbit-edge backend should exploit that
  sparsity.

  Why I removed the 1D specializer

  The active-shift specializer did exactly this, but only for 1D periodic
  systems. You correctly pointed out that this targets the example too closely.
  The proper version must be generic over lattice dimension and periodic flags.

  Correct next theoretical/implementation target

  You are right that the commutator must become more TS-native in cost.

  The generic target is:

  orbit_edges(Ha, q)

  implemented with dimension-independent active relative shifts, not full
  all_shifts.

  Then:

  orbit_liouvillian(Ha, O)

  should apply that edge generator directly without making per-q temporary
  operators.

  In other words:

  current:
      orbit_liouvillian -> commutator -> all relative shifts

  desired:
      orbit_liouvillian -> active orbit-edge generator -> only contributing
  relative shifts

  And it must be generic, not MFIM/1D-specific.

  Bottom line

  You are correct in spirit:

  ▎ If the orbit Liouvillian still scans all translations, then the method is 
  ▎ not fully exploiting TS structure.

  The current method is TS-native in representation and splitting, but the
  Liouvillian backend is still too close to “resum over relative shifts.” The
  next real work is a generic active-shift Pauli-orbit edge generator.

✻ Brewed for 36s

❯ Yes i agree. but before you conitnue, make a markdown flow tree of all the    
  things we have attempted, and why they are failing. This for internal use, so 
  we can verify we aren't going in circles                                      

● Write(benchmark/ts_trotter_attempts_flow.md)
  ⎿  Wrote 159 lines to benchmark/ts_trotter_attempts_flow.md
       1 # TS Trotter attempts flow
       2 
       3 Internal tracking notes to avoid going in circles.
       4 
       5 ## Root problem
       6 
       7 Current TS Trotter does:
       8 
       9 ```text
      10 H::OperatorTS -> resum(H) -> ordinary translated gates
     … +192 lines (ctrl+o to expand)

● Wrote the internal flow-tree document:

  benchmark/ts_trotter_attempts_flow.md

  It tracks each attempted avenue, what happened, and why it failed or remains
  useful. It ends with the current understanding and the real next target: a
  generic dimension-independent active relative-shift generator for Pauli-orbit
  edges.

✻ Churned for 51s

※ recap: We’re building a generic TS-native second-order Trotter method using orbit Liouvillian flows rather than expanded translated gates. Next step is to design a dimension-independent active relative-shift generator and validate it against the current commutator.

❯ good. now lets attemp to find a generic active-shift pauli orbit edge         
  generator                                                                     
  ⎿  Interrupted · What should Claude do instead?

❯ good. now lets attemp to find a generic active-shift pauli orbit edge         
  generator                                                                     
  ⎿  Interrupted · What should Claude do instead?

❯ good. now lets attemp to find a generic active-shift pauli orbit edge         
  generator                                                                     
● I’ll implement a generic active-shift edge generator for orbit_edges:
  dimension-independent, using support bit iteration and a bitset for
  deduplicating relative shifts, with no 1D-specific assumptions. Then I’ll
  validate exactly against the existing commutator and measure edge-generation
  time/allocations.

  Read 1 file (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Added 66 lines, removed 5 lines
      257      return op
      258  end
      259  
      260 +@inline function _coord(site::Integer, Ls::Tuple, dim::Integer
          +)                                                             
      261 +    x = site - 1                                              
      262 +    for k in 1:(dim - 1)                                      
      263 +        x = div(x, Ls[k])                                     
      264 +    end                                                       
      265 +    return mod(x, Ls[dim]) + 1                                
      266 +end                                                           
      267 +                                                              
      268 +@inline function _linear_site(coords::Tuple, Ls::Tuple)       
      269 +    site = 1                                                  
      270 +    stride = 1                                                
      271 +    for k in eachindex(Ls)                                    
      272 +        site += (coords[k] - 1) * stride                      
      273 +        stride *= Ls[k]                                       
      274 +    end                                                       
      275 +    return site                                               
      276 +end                                                           
      277 +                                                              
      278 +function _active_shift_tuple(site_p::Integer, site_q::Integer,
          + Ls::Tuple, Ps::Tuple)                                        
      279 +    valid = true                                              
      280 +    shifts = ntuple(length(Ls)) do k                          
      281 +        cp = _coord(site_p, Ls, k)                            
      282 +        cq = _coord(site_q, Ls, k)                            
      283 +        if Ps[k]                                              
      284 +            d = mod(cp - cq, Ls[k])                           
      285 +            iszero(d) ? Ls[k] : d                             
      286 +        else                                                  
      287 +            valid &= cp == cq                                 
      288 +            1                                                 
      289 +        end                                                   
      290 +    end                                                       
      291 +    return valid, shifts                                      
      292 +end                                                           
      293 +                                                              
      294 +function _active_shift_key(shifts::Tuple, Ls::Tuple)          
      295 +    return _linear_site(shifts, Ls)                           
      296 +end                                                           
      297 +                                                              
      298  function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStr
           ingTS; epsilon::Real=0, maxlength::Int=1000)
      299      checklength(A, q)
      300      Ls = qubitsize(A)
      301      Ps = periodicflags(A)
      302 +    N = Base.prod(Ls)                                         
      303      d = emptydict(A)
      304      qrep = representative(q)
      305 +    uq = qrep.v | qrep.w                                      
      306 +    iszero(uq) && return typeof(A)(collect(keys(d)), collect(v
          +alues(d)))                                                    
      307 +                                                              
      308      for (p, c) in zip(A.strings, A.coeffs)
      309          prep = representative(p)
      268 -        for s in all_shifts(paulistringtype(A))               
      269 -            out, k = commutator(prep, shift(qrep, Ls, Ps, s)) 
      270 -            coeff = c * k                                     
      271 -            if (k != 0) && (abs(coeff) > epsilon) && pauli_wei
          -ght(out) < maxlength                                          
      272 -                setwith!(+, d, PauliStringTS{Ls,Ps}(out), coef
          -f)                                                            
      310 +        up = prep.v | prep.w                                  
      311 +        iszero(up) && continue                                
      312 +        seen = falses(N)                                      
      313 +        xp = up                                               
      314 +        while !iszero(xp)                                     
      315 +            site_p = trailing_zeros(xp) + 1                   
      316 +            xp &= xp - one(xp)                                
      317 +                                                              
      318 +            xq = uq                                           
      319 +            while !iszero(xq)                                 
      320 +                site_q = trailing_zeros(xq) + 1               
      321 +                xq &= xq - one(xq)                            
      322 +                                                              
      323 +                valid, shifts = _active_shift_tuple(site_p, si
          +te_q, Ls, Ps)                                                 
      324 +                valid || continue                             
      325 +                key = _active_shift_key(shifts, Ls)           
      326 +                seen[key] && continue                         
      327 +                seen[key] = true                              
      328 +                                                              
      329 +                out, k = commutator(prep, shift(qrep, Ls, Ps, 
          +shifts))                                                      
      330 +                coeff = c * k                                 
      331 +                if (k != 0) && (abs(coeff) > epsilon) && pauli
          +_weight(out) < maxlength                                      
      332 +                    setwith!(+, d, PauliStringTS{Ls,Ps}(out), 
          +coeff)                                                        
      333 +                end                                           
      334              end
      335          end
      336      end

  Ran 1 shell command (ctrl+o to expand)

● Implemented a first generic active-shift orbit_edges.

  What changed

  orbit_edges(A::OperatorTS, q::PauliStringTS) no longer loops over all
  translations.

  Instead, for each Hamiltonian representative support site and each
  input-string support site, it computes the unique relative shift that aligns
  them:

  site_p in support(P)
  site_q in support(Q)
    -> relative shift δ such that shift(Q, δ) overlaps P at site_p

  It deduplicates relative shifts with a Boolean vector of length prod(Ls) and
  evaluates only those candidate shifts.

  This is generic over:

  - arbitrary Ls dimensionality,
  - periodic and non-periodic Ps,
  - no MFIM assumptions,
  - no 1D-only branch.

  Validation

  Checked against the old mathematical reference:

  orbit_edges(Hts, q) == commutator(Hts, OperatorTS(q))

  via resummed norm for:

  N=6, Ls=(6,)
  N=8, Ls=(8,)
  N=6, Ls=(2,3)

  All passed:

  (N = 6, Ls = (6,), ok = true)
  (N = 8, Ls = (8,), ok = true)
  (N = 6, Ls = (2, 3), ok = true)

  Caveat

  This currently uses:

  seen = falses(N)

  inside orbit_edges, so it allocates one Boolean vector per Hamiltonian term
  per input orbit. It is generic and correct, but not yet allocation-optimized.

  Next step is to measure:

  - candidate shift count reduction,
  - orbit_edges time/allocations,
  - whether using this active edge generator in a future direct graph backend is
  promising.

✻ Brewed for 1m 59s

❯ continue                                                                      

  Ran 1 shell command (ctrl+o to expand)

● Measured the generic active-shift orbit_edges.

  Edge counts

  The folded edge counts are unchanged, as expected, because the new
  implementation is mathematically equivalent.

  N=32, M=12, lenO=4096

  ALL folded_edges = 45569
  avg per input/orbit = 3.708
  max = 8

  N=64, M=14, lenO=16384

  ALL folded_edges = 208410
  avg per input/orbit = 4.240
  max = 9

  Allocation in orbit_edges

  Compared to the old generic full-shift orbit_edges, allocations are mixed:

  N=32

  orbit 1 alloc = 13.4 MB
  orbit 2 alloc = 4.46 MB
  orbit 3 alloc = 4.46 MB

  N=64

  orbit 1 alloc = 24.9 MB
  orbit 2 alloc = 21.3 MB
  orbit 3 alloc = 20.9 MB

  The per-input seen = falses(N) allocation is likely showing up here,
  especially for larger N.

  Shift sparsity confirmed

  The active-shift idea is theoretically justified:

  N=32

  full shifts = 393216
  overlap candidates = 76775  = 19.52%
  nonzero commutators = 45569 = 11.59%

  N=64

  full shifts = 3145728
  overlap candidates = 350724 = 11.15%
  nonzero commutators = 208410 = 6.63%

  So active-shift enumeration reduces candidate work by roughly:

  5x at N=32
  9x at N=64

  before folding.

  Current issue

  The generic active-shift implementation is correct and dimension-independent,
  but still allocates in hot usage because it creates a fresh Boolean seen
  vector per input orbit.

  Next improvement should be allocation-focused but generic:

  replace seen = falses(N)

  with a reusable integer stamp array or workspace object, likely internal to a
  future orbit_edges! / orbit_liouvillian! style API.

  That would keep the generic active-shift logic while avoiding per-input
  allocations.

✻ Brewed for 52s

※ recap: We’re building a generic TS-native Trotter method using orbit-level Liouvillian exponentials instead of expanded translated gates. Next, make the generic active-shift orbit edge generator allocation-light by replacing per-call Boolean arrays with reusable workspace/stamping.

❯ before that, benchmark against rk4                                            

● Write(~/.claude/jobs/513a65b3/tmp/compare_active_edges_rk4.jl)
  ⎿  Wrote 32 lines to
     ../../../.claude/jobs/513a65b3/tmp/compare_active_edges_rk4.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function bench(N,M,steps; reps=3)
     … +30 lines (ctrl+o to expand)

  Ran 4 shell commands (ctrl+o to expand)

● Benchmarked the generic active-shift orbit_edges state against RK4/expanded
  Trotter, using current production TrotterTS(krylovdim=4).

  Important: TrotterTS hot path still uses orbit_liouvillian = commutator, so
  active orbit_edges does not accelerate TrotterTS yet. It only affects graph
  inspection/future backends.

  Runtime benchmark

  MFIM, dt=0.1, no fout, best of 3.

  ┌───────────────┬─────────────┬─────────────┬────────────────┬──────────┐
  │     Case      │     RK4     │  Expanded   │ TrotterTS(k=4) │  TS/RK4  │
  │               │             │   Trotter   │                │          │
  ├───────────────┼─────────────┼─────────────┼────────────────┼──────────┤
  │ N=16, M=10,   │   0.0536s,  │   0.0685s,  │       0.327s,  │   6.11x  │
  │ steps=20      │     0.063GB │     0.116GB │        0.692GB │   slower │
  ├───────────────┼─────────────┼─────────────┼────────────────┼──────────┤
  │ N=32, M=12,   │    0.319s,  │    0.522s,  │       1.652s,  │   5.17x  │
  │ steps=20      │     0.270GB │     0.900GB │        2.916GB │   slower │
  ├───────────────┼─────────────┼─────────────┼────────────────┼──────────┤
  │ N=64, M=14,   │    0.775s,  │    1.839s,  │       0.250s,  │    0.32x │
  │ steps=10      │     0.478GB │     3.577GB │        0.424GB │          │
  └───────────────┴─────────────┴─────────────┴────────────────┴──────────┘

  The N=64 timing looks surprisingly favorable because final TrotterTS length
  was much smaller:

  RK4 final len:       16384
  TrotterTS final len: 5888

  So I also checked observable histories.

  Observable comparison vs RK4

  Using normalized trace_product(O0, O; scale=1):

  ┌─────────────────┬──────────────┬───────────┬──────────┬──────────────┐
  │      Case       │ rel history  │ max abs   │  final   │    final     │
  │                 │     diff     │   diff    │   RK4    │  TrotterTS   │
  ├─────────────────┼──────────────┼───────────┼──────────┼──────────────┤
  │ N=32, M=12,     │      0.00292 │    0.0762 │  11.6627 │      11.6124 │
  │ steps=20        │              │           │          │              │
  ├─────────────────┼──────────────┼───────────┼──────────┼──────────────┤
  │ N=64, M=14,     │      0.00274 │    0.1548 │  17.0094 │      16.9476 │
  │ steps=10        │              │           │          │              │
  └─────────────────┴──────────────┴───────────┴──────────┴──────────────┘

  So TrotterTS is not wildly wrong in this observable, but the shorter final
  representation means we should be cautious comparing wall time directly. It
  may be effectively more dissipative/truncating in Pauli space.

  Takeaway

  - Current generic TrotterTS(krylovdim=4) is usually slower than RK4 at modest
  sizes.
  - It can look faster when it produces much smaller truncated operators, so
  final runtime alone is not enough.
  - Active orbit_edges is correct and shows graph sparsity, but it is not yet
  used in the hot Liouvillian path because the direct orbit_edges backend
  allocated too much.
  - The real missing piece remains: a generic active-shift orbit_liouvillian
  backend that avoids both:
    a. scanning all shifts like the current commutator,
    b. creating per-q temporary OperatorTS edge objects like the naive graph
  backend.

✻ Churned for 3m 29s

❯ For now, use orbit_edges in the hot path. I want to verify that we are        
  actually doing less work compared to RK4. I understand that orbit splitting   
  reduces from Nm to m, but if the orbit update is order N slower than direct   
  updates, then we don't gain anything                                          

● I’ll route orbit_liouvillian back through orbit_edges in the hot path, then
  benchmark/instrument against RK4 so we can compare actual work: orbit-edge
  calls, shifted candidates, allocations, and wall time.

● Update(src/translation_symmetry.jl)
  ⎿  Added 10 lines, removed 1 line
      339  
      340  function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
           gth::Int=1000)
      341      checklength(A, B)
      342 -    return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
          -axlength) / hbar                                              
      342 +    d = emptydict(A)                                          
      343 +    for (q, cq) in zip(B.strings, B.coeffs)                   
      344 +        edges = orbit_edges(A, q; epsilon=epsilon, maxlength=m
          +axlength)                                                     
      345 +        for (p, c) in zip(edges.strings, edges.coeffs)        
      346 +            setwith!(+, d, p, cq * c)                         
      347 +        end                                                   
      348 +    end                                                       
      349 +    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      350 +    o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) 
          +: o                                                           
      351 +    return 1im * o / hbar                                     
      352  end
      353  
      354  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)

● Write(~/.claude/jobs/513a65b3/tmp/edge_hotpath_work_compare.jl)
  ⎿  Wrote 56 lines to
     ../../../.claude/jobs/513a65b3/tmp/edge_hotpath_work_compare.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +62 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Routed orbit_liouvillian through orbit_edges in the hot path and benchmarked
  against RK4.

  Correctness

  test/evolution.jl still passes:

  evolve | 8 passed

  Work comparison on same input state

  This uses an evolved state as the input, then measures one step.

  N=16, M=10

  start_len = 1024
  full relative-shift checks = 49,152
  folded orbit_edges = 9,742
  avg edges/input/orbit = 3.17
  edge_count_time = 0.00094s
  edge_count_alloc = 0.013 GB

  Runtime:

  ┌───────────┬──────────┬──────────┬───────────┐
  │  Method   │   Time   │  Alloc   │ Final len │
  ├───────────┼──────────┼──────────┼───────────┤
  │ RK4       │ 0.00304s │ 0.108 GB │      1024 │
  ├───────────┼──────────┼──────────┼───────────┤
  │ TrotterTS │  0.0214s │ 0.385 GB │      1024 │
  └───────────┴──────────┴──────────┴───────────┘

  N=32, M=12

  start_len = 4096
  full relative-shift checks = 393,216
  folded orbit_edges = 45,569
  avg edges/input/orbit = 3.71
  edge_count_time = 0.00494s
  edge_count_alloc = 0.017 GB

  Runtime:

  ┌───────────┬─────────┬──────────┬───────────┐
  │  Method   │  Time   │  Alloc   │ Final len │
  ├───────────┼─────────┼──────────┼───────────┤
  │ RK4       │ 0.0181s │ 0.020 GB │      4096 │
  ├───────────┼─────────┼──────────┼───────────┤
  │ TrotterTS │  0.127s │ 0.300 GB │      4096 │
  └───────────┴─────────┴──────────┴───────────┘

  N=64, M=14

  start_len = 16384
  full relative-shift checks = 3,145,728
  folded orbit_edges = 208,410
  avg edges/input/orbit = 4.24
  edge_count_time = 0.0308s
  edge_count_alloc = 0.063 GB

  Runtime:

  ┌───────────┬────────┬──────────┬───────────┐
  │  Method   │  Time  │  Alloc   │ Final len │
  ├───────────┼────────┼──────────┼───────────┤
  │ RK4       │ 0.129s │ 0.088 GB │     16384 │
  ├───────────┼────────┼──────────┼───────────┤
  │ TrotterTS │ 0.252s │ 0.504 GB │     16384 │
  └───────────┴────────┴──────────┴───────────┘

  Interpretation

  The orbit-edge graph is doing less combinatorial work:

  N=64:
  full shift checks: 3.15M
  folded edges:      208k
  avg degree:        ~4.24

  So the graph is sparse and the theoretical reduction is real.

  But the current TrotterTS hot path is still slower because Krylov uses many
  orbit-edge/Liouvillian applications and creates many intermediate OperatorTS
  objects.

  For m=3, krylovdim=4, one TrotterTS step does roughly:

  5 orbit flows * 4 Krylov applications = 20 orbit Liouvillian applications

  RK4 does:

  4 full-H commutators

  So even if each orbit-edge application is sparse, the Krylov machinery
  multiplies the work.

  Important conclusion

  Using orbit_edges in the hot path verifies we are reducing graph degree
  relative to scanning all shifts, but this alone is not enough to beat RK4
  because the current orbit-flow exponential backend is too expensive.

  The problem is now clearly:

  orbit exponential action implementation

  not:

  orbit graph sparsity

  The graph is sparse; the Krylov wrapper is expensive.

✻ Cogitated for 1m 40s

❯ Then why are we not focussing on that?                                        

● You’re right. We should focus there now.

  The mistake in the last several loops was treating orbit_edges / Liouvillian
  application as the main missing piece. It is important, but the measurements
  now say:

  orbit graph sparsity: good
  orbit_edges generation: not the dominant cost
  Krylov orbit-flow wrapper: dominant cost

  So the next target should be the orbit exponential action, not the edge
  generator.

  What the data says

  For one TrotterTS(krylovdim=4) step with m=3:

  orbit factors = 5
  Krylov dim = 4
  Liouvillian applications ≈ 20
  inner products ≈ 100
  residual updates ≈ 50
  truncations ≈ 75

  RK4 does:

  4 full-H commutators

  Even if each orbit-edge application is sparse, Krylov wraps it in a lot of
  extra machinery.

  So the current method is mathematically clean but computationally overbuilt
  for the target step sizes.

  The actual question

  We need a cheaper way to approximate/apply:

  e^{τ L_a} O

  than generic Arnoldi Krylov.

  That means we need to focus on orbit-flow backends, e.g.

  1. Low-dimensional local component exponentials, but implemented as reusable
  plans, not rebuilt every step.
  2. Lanczos/short recurrence if the orbit Liouvillian is skew-Hermitian under 
  the right inner product, reducing orthogonalization from O(k²) to O(k).
  3. Fixed low-order symmetric exponential approximants, e.g. Cayley/Padé, if
  they preserve enough structure and are cheaper than Krylov.
  4. A direct Pauli-orbit propagation formula on small graph components.

  Why Krylov is expensive

  The current Arnoldi version pays:

  for each orbit flow:
    build basis V1...Vk as OperatorTS
    compute many trace_product inner products
    orthogonalize each new vector against all previous ones
    repeatedly construct residual OperatorTS objects

  For k=4, the overhead is already comparable to or larger than the Liouvillian
  application.

  So making orbit_edges faster will not solve the fundamental issue unless the
  orbit-flow backend also stops multiplying the number of Liouvillian
  applications and operator algebra calls.

  What we should do next

  I would pivot to a planned local-component backend, but do it correctly:

  Do not do this again

  for every step:
    discover component
    build matrix
    exp(matrix)

  That failed.

  Instead

  Build and cache orbit-flow plans:

  struct OrbitComponentPlan
      basis::Vector{PauliStringTS}
      index::Dict{PauliStringTS,Int}
      A::Matrix{ComplexF64}
      expA_dt::Matrix{ComplexF64}
  end

  For each Hamiltonian orbit Ha and time step dt, cache components as they are
  discovered.

  Then an orbit flow is:

  for each active q:
      find component plan
      gather coeffs
      coeffs2 = expA_dt * coeffs
      scatter back

  This avoids repeated graph discovery and repeated dense exponentiation.

  Why this is the right focus

  The component measurements were promising:

  N=32:
  p50 component size: 4–8
  p90 component size: 16–32
  max mostly ≤ 128

  These are exactly sizes where local dense exponentials can be practical if 
  cached.

  This would make orbit flow closer to:

  one cached matrix-vector multiply per component

  instead of:

  20 Liouvillian applications + 100 inner products per step

  So yes: next step

  We should implement/prototype cached component plans for one orbit flow and
  measure:

  plan build time
  number of components
  component sizes
  exp build time
  per-step apply time
  allocations
  correctness vs Krylov/dense

  If cached component apply is fast, then it becomes the real TS-native orbit
  exponential backend.

✻ Brewed for 38s

❯ Ok so we are changing the krylov backend, because theoretically the krylov    
  update does too much work, even if the component actions are O(1) operations? 

● Yes — exactly.

  More precisely:

  The orbit graph action may be cheap, but Krylov asks for it too many times

  For one orbit flow:

  O \mapsto e^{τ L_a} O,

  Krylov does not apply this exponential directly. It builds a Krylov basis:

  O,\; L_a O,\; L_a^2 O,\; L_a^3 O,\ldots

  So even if one L_a O application is cheap, Krylov with dimension k needs
  about:

  k orbit-Liouvillian applications
  +
  O(k²) inner products
  +
  O(k²) residual updates
  +
  many truncations

  Then TrotterTS applies this for every Strang orbit factor.

  For m = 3, second-order Strang has:

  2m - 1 = 5 orbit flows

  With krylovdim = 4:

  5 * 4 = 20 orbit-Liouvillian applications per Trotter step

  plus roughly:

  100 trace_product inner products
  50 residual updates
  75 truncations

  So even if the orbit graph action were O(1) or very sparse, the Krylov wrapper
  multiplies it heavily.

  Why component actions are different

  If the orbit Liouvillian graph decomposes into small connected components,
  then for each component we can precompute:

  E_a = e^{τ A_a}

  where A_a is the finite matrix representation of L_a on that component.

  Then each orbit flow is just:

  gather coefficients
  matrix-vector multiply with cached exp(A)
  scatter coefficients

  No repeated Krylov basis construction. No repeated orthogonalization. No k
  applications of L_a.

  So if component sizes are small and reusable, component exponentials can turn
  one orbit flow from:

  k Liouvillian applications + Krylov overhead

  into:

  one cached local matrix-vector multiply

  That is the theoretical advantage we are chasing.

  The key condition

  This only works if the orbit graph components are actually small or can be
  capped/truncated sensibly.

  Our sampling suggested:

  p50 component size ~ 4–8
  p90 component size ~ 16–32

  which is promising.

  But the naive component implementation failed because it rebuilt components
  and exponentials every time.

  Therefore

  Yes, we are changing the Krylov backend because:

  Krylov is a generic exponential-action method,
  but it does too much repeated work for this sparse orbit graph.

  The more specialized TS-native backend should exploit the orbit graph’s
  connected components directly.

✻ Baked for 25s

❯ why do we need krylov in the first place? if we are attempting to make a TS   
  Trotter, just a TS commutator should be enough right?                         

● A single TS commutator is enough to compute the derivative, not the Trotter 
  substep.

  That distinction is the core issue.

  What a TS commutator gives

  For one orbit Hamiltonian:

  H_a^{orb}

  the TS commutator gives:

  L_a(O) = \frac{i}{\hbar}[H_a^{orb}, O].

  This is the infinitesimal generator:

  \frac{dO}{dt} = L_a(O).

  So yes, a TS commutator is enough to know the instantaneous direction of
  motion.

  But a Trotter factor is not the derivative. It is the finite-time flow:

  O(t+\tau) = e^{\tau L_a} O(t).

  To apply a Trotter step, we need the exponential action of the orbit
  Liouvillian.

  Non-TS Trotter gets this exponential for free

  For a single Pauli generator P, the finite-time action is closed-form:

  e^{i\theta P/2} Q e^{-i\theta P/2}
  =
  \begin{cases}
  Q, & [P,Q]=0, \\\\
  \cos(\theta) Q + \sin(\theta) Q', & \{P,Q\}=0.
  \end{cases}

  That is why non-TS Trotter is fast: one Pauli rotation is an exact finite-time
  exponential action.

  TS orbit generator is not a single Pauli

  For a TS orbit:

  H_a^{orb} = \sum_r T_r P_a T_r^\dagger,

  the generator is a sum of many translated Pauli strings.

  Even though the commutator is TS-native:

  L_a(O) = i[H_a^{orb}, O],

  the finite flow is:

  e^{\tau L_a} O
  =
  O + \tau L_a O
  + \frac{\tau^2}{2}L_a^2O
  + \frac{\tau^3}{6}L_a^3O
  + \cdots.

  A single commutator gives only the first term.

  If we use one commutator, what method is that?

  If an orbit substep were:

  O \leftarrow O + \tau L_a O,

  that is explicit Euler for the orbit flow.

  If we add:

  O + \tau L_a O + \frac{\tau^2}{2}L_a^2O,

  that is Taylor-2.

  Those are commutator-based integrators, not exact Trotter orbit exponentials.
  Taylor-2 was the thing we rejected because it becomes RK-like rather than a
  true Trotter primitive.

  Why Krylov appeared

  Krylov was one generic way to approximate:

  e^{\tau L_a} O

  without forming the full matrix of L_a.

  It was conceptually correct:

  orbit Trotter factor = exponential of orbit Liouvillian

  but too expensive in the current implementation.

  What we actually want instead of Krylov

  We want a TS analogue of the single-Pauli rotation rule.

  For one orbit generator, the orbit Liouvillian defines a graph:

  [Q] -> [R1], [R2], ...

  The finite orbit flow is the exponential of this graph’s adjacency/Liouvillian
  matrix.

  If the connected graph component of [Q] is small, we can exactly exponentiate
  that small component:

  e^{\tau A_{\mathrm{component}}}

  That would be the TS-orbit analogue of the two-dimensional single-Pauli
  rotation.

  So the hierarchy is:

  single commutator       = derivative only
  Krylov                  = generic exponential action
  component exponential   = specialized TS Pauli-orbit exponential

  Bottom line

  A TS commutator is enough for RK-style time stepping.

  A TS Trotter factor needs a finite-time orbit exponential:

  e^{\tau L_a}.

  Krylov was a generic placeholder for that exponential. The better target is a
  small-component Pauli-orbit exponential that uses the TS graph structure
  directly.

✻ Baked for 30s

❯ ok then lets directly focus on finding a TS pauli rotation rule               
● Yes. That is the right target.

  We need the TS analogue of the single-Pauli rotation rule.

  For an ordinary Pauli generator P, the key fact is:

  \operatorname{ad}_P^2(Q) \propto Q,

  so the orbit of Q under repeated commutators closes in a 2D space:

  \mathrm{span}\{Q,\; PQ\}.

  That gives the closed-form rotation:

  Q \mapsto \cos(\theta)Q + \sin(\theta)Q'.

  For a TS orbit generator:

  [P] = \sum_r T_r P T_r^\dagger,

  we need to understand the closure of:

  \operatorname{ad}_{[P]}([Q]).

  The candidate TS Pauli rotation rule is:

  ▎ For each input orbit [Q], build the connected component of the Pauli-orbit 
  ▎ graph generated by [P]; the exact orbit rotation is the matrix exponential 
  ▎ on that component.

  So the “rule” is no longer generally 2D. It is:

  [Q] \mapsto \sum_{R \in C(Q)} \left(e^{\theta A_C}\right)_{R,Q} [R],

  where C(Q) is the connected component under the orbit-edge graph of [P].

  This is the generic TS Pauli rotation rule.

  Derivation

  Let:

  [P] = \sum_r P_r,
  \qquad
  [Q] = \sum_x Q_x.

  The orbit Liouvillian is:

  L_P([Q]) = i[[P],[Q]].

  Using translation covariance, this folds to relative shifts:

  L_P([Q])
  =
  i \sum_\delta [P_\delta, Q]_{\mathrm{folded}}.

  For each relative shift δ:

  - if P_δ commutes with Q, no edge;
  - if not, the commutator gives a Pauli string:

  R = P_\delta Q

  up to phase;
  - fold R to its canonical TS orbit representative [R];
  - accumulate the phase/multiplicity.

  Thus:

  L_P [Q] = \sum_R A_{R,Q} [R].

  Now repeated commutators are just repeated graph transitions:

  L_P^2[Q], L_P^3[Q], \ldots

  So the smallest invariant subspace containing [Q] is the connected component
  of Q in this graph.

  Then the exact TS orbit rotation is:

  e^{\tau L_P}[Q]
  =
  \sum_{R \in C(Q)}
  \left(e^{\tau A_C}\right)_{R,Q}
  [R].

  This is the TS version of Pauli rotation.

  Relation to ordinary Pauli rotation

  For a single ordinary Pauli P, the component size is at most 2:

  Q <-> PQ

  So:

  A_C

  is a 2×2 skew-Hermitian matrix, and the exponential reduces to sine/cosine.

  For a TS orbit generator [P], the component may have size:

  2, 4, 8, 16, ...

  depending on overlaps and folding. The rotation is still exact, but the
  component matrix may be larger than 2×2.

  So the generalization is:

  ordinary Pauli rotation:
      exact 2-node component exponential

  TS Pauli-orbit rotation:
      exact connected-component exponential

  What we need to verify

  The key question is whether components remain small enough in practice.

  We already sampled and found promising sizes:

  p50 ~ 4–8
  p90 ~ 16–32

  But now we need to make this rigorous enough for implementation:

  1. For a given orbit generator [P], define edge matrix A.
  2. For each active [Q], compute full connected component C(Q).
  3. Verify closure:

  L_P(C) \subseteq C.
  4. Build A_C.
  5. Apply:

  e^{\tau A_C}
  6. Compare against dense exact orbit evolution for small N.

  Implementation shape

  This should become the orbit flow primitive:

  pauli_orbit_rotation(Ha, O, τ)

  Internally:

  for each unassigned active string q in O:
      C = connected_component(Ha, q)
      A = orbit_component_matrix(Ha, C)
      coeffs = gather O coefficients on C
      coeffs′ = exp(τA) * coeffs
      scatter coeffs′

  Important: this is the conceptual implementation. The first correct version
  can rebuild components. The optimized version caches (Ha, component) matrices.

  Why this is not RK4

  Because this applies:

  e^{\tau L_P}

  exactly on each orbit component.

  It is not a Taylor polynomial in L_P. It is the exact finite-time orbit
  rotation, just like ordinary Trotter applies exact finite-time single-Pauli
  rotations.

  Why this is not MFIM-specific

  This uses only:

  translation orbit structure
  Pauli commutation graph
  canonical TS folding

  It does not assume commuting translated terms, 1D geometry, or even/odd
  coloring.

  Next step

  I would implement a clean internal function for the exact component rotation,
  not optimized:

  _orbit_rotation_component(Ha, O, τ; componentmax=...)

  Then test it against dense exact for a single orbit flow:

  Ha = one TS orbit term
  O  = TS operator
  τ  = small time

  Matrix exact:
  U = exp(im * Matrix(resum(Ha)) * τ)
  O_exact = U * Matrix(O) * U'

  Component:
  O_comp = _orbit_rotation_component(Ha, O, τ)

  compare Matrix(O_comp) to O_exact

  If this passes, we have found the TS Pauli rotation rule. Performance
  optimization comes after.

✻ Brewed for 39s

❯ /effort high                                                                  
  ⎿  Set effort level to high: Comprehensive implementation with extensive
     testing and documentation

❯ ok try                                                                        

  Read 2 files, ran 1 shell command (ctrl+o to expand)

● I’ll try the exact Pauli-orbit rotation rule in a scratch prototype first:
  build the closed orbit graph component for a single orbit Hamiltonian,
  exponentiate the finite Liouvillian matrix, and compare against dense exact
  evolution under resum(Ha). If this works, we can decide how to integrate it
  cleanly.

● Write(~/.claude/jobs/513a65b3/tmp/orbit_rotation_exact_test.jl)
  ⎿  Wrote 106 lines to
     ../../../.claude/jobs/513a65b3/tmp/orbit_rotation_exact_test.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function TFIM(N,h)
       6     H=Operator(N)
       7     H += -h,"X",1
       8     H += "Z",1,"Z",2
       9     OperatorTS{(N,)}(H)
      10 end
     … +101 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/orbit_rotation_exact_test.jl)
  ⎿  Added 1 line, removed 1 line
      75          coeff[index[q]] += c
      76      end
      77      coeff2=exp(τ*A)*coeff
      78 -    return typeof(O)(basis,coeff2), length(basis), nnz(A)      
      78 +    return typeof(O)(basis,coeff2), length(basis), count(!iszer
         +o, A)                                                          
      79  end
      80  
      81  function dense_orbit_rotate(Ha,O,τ; hbar=1)

  Ran 1 shell command (ctrl+o to expand)

● I tested the exact TS Pauli-orbit rotation rule in a scratch prototype.

  Prototype

  For a single orbit Hamiltonian Ha and an input O, the prototype:

  1. builds the closure of the TS Pauli-orbit graph using orbit_edges(Ha, q),
  2. constructs the finite matrix A for the orbit Liouvillian on that closed
  component,
  3. applies:

  exp(τ * A) * coeffs
  4. compares the resulting OperatorTS against dense exact evolution under:

  resum(Ha)

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/orbit_rotation_exact_test.jl

  Results

  TFIM, N=4

  orbit 1, τ=0.2: component_size=12, nnzA=24, rel_dense_err=1.66e-16
  orbit 1, τ=0.1: component_size=12, nnzA=24, rel_dense_err=1.82e-16

  orbit 2, τ=0.2: component_size=6, nnzA=6, rel_dense_err=1.66e-16
  orbit 2, τ=0.1: component_size=6, nnzA=6, rel_dense_err=2.31e-16

  TFIM, N=6

  orbit 1, τ=0.1:  component_size=12, nnzA=24, rel_dense_err=1.98e-16
  orbit 1, τ=0.05: component_size=12, nnzA=24, rel_dense_err=1.34e-16

  orbit 2, τ=0.1:  component_size=7, nnzA=10, rel_dense_err=2.93e-16
  orbit 2, τ=0.05: component_size=7, nnzA=10, rel_dense_err=2.63e-16

  MFIM, N=6

  orbit 1, τ=0.1: component_size=12, nnzA=24, rel_dense_err=1.98e-16
  orbit 2, τ=0.1: component_size=7,  nnzA=10, rel_dense_err=2.93e-16
  orbit 3, τ=0.1: component_size=6,  nnzA=6,  rel_dense_err=1.10e-16

  Conclusion

  This verifies the mathematical TS Pauli-orbit rotation rule:

  exact orbit flow = dense exponential on the closed Pauli-orbit graph component

  For small systems, it matches dense exact evolution to machine precision.

  This is the real TS analogue of ordinary single-Pauli rotation:

  ordinary Pauli rotation:
    component size ≤ 2

  TS Pauli-orbit rotation:
    finite connected component in Pauli-orbit graph

  The next implementation step is to turn this prototype into an internal
  orbit-flow backend, with caching/planning so we do not rebuild the component
  and dense matrix every step.

✻ Cooked for 1m 58s

❯ yes go ahead                                                                  

● Update(src/evolution.jl)
  ⎿  Added 7 lines, removed 5 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, krylovdim=4)                          
      29 +    TrotterTS(; order=2, krylovdim=4, componentmax=1024)       
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34 -Liouvillians. Each orbit flow is applied as a matrix-free Arnol
         -di/Krylov exponential                                          
      35 -action with dimension `krylovdim`.                             
      34 +Liouvillians. Each orbit flow is applied by exponentiating conn
         +ected components of                                            
      35 +the Pauli-orbit graph. Components larger than `componentmax` fa
         +ll back to a                                                   
      36 +matrix-free Arnoldi/Krylov exponential action with dimension `k
         +rylovdim`.                                                     
      37  """
      38  struct TrotterTS <: AbstractEvolutionMethod
      39      order::Int
      40      krylovdim::Int
      41 +    componentmax::Int                                          
      42  end
      41 -TrotterTS(; order::Integer=2, krylovdim::Integer=4) =          
      42 -    TrotterTS(Int(order), Int(krylovdim))                      
      43 +TrotterTS(; order::Integer=2, krylovdim::Integer=4, componentma
         +x::Integer=1024) =                                             
      44 +    TrotterTS(Int(order), Int(krylovdim), Int(componentmax))   
      45  
      46  """
      47      RK4()

● I’ll integrate the exact component-exponential orbit flow as the default
  TrotterTS backend, with Krylov fallback only if a component exceeds a cap.
  Then I’ll test against dense exact on small systems and benchmark.

● Update(src/evolution.jl)
  ⎿  Added 72 lines, removed 1 line
      347      return truncation(_axpy_combination(V, coeffs, O))
      348  end
      349  
      350 +function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
          +:PauliStringTS, componentmax::Int, maxlength::Int)            
      351 +    component = typeof(seed)[seed]                            
      352 +    index = Dict{typeof(seed),Int}(seed => 1)                 
      353 +    queue_index = 1                                           
      354 +    while queue_index <= length(component)                    
      355 +        q = component[queue_index]                            
      356 +        queue_index += 1                                      
      357 +        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      358 +        for r in edges.strings                                
      359 +            haskey(index, r) && continue                      
      360 +            length(component) >= componentmax && return compon
          +ent, index, true                                              
      361 +            push!(component, r)                               
      362 +            index[r] = length(component)                      
      363 +        end                                                   
      364 +    end                                                       
      365 +    return component, index, false                            
      366 +end                                                           
      367 +                                                              
      368 +function _orbit_component_matrix(Ha::Operator{<:PauliStringTS}
          +, component, index, hbar::Real, maxlength::Int)               
      369 +    n = length(component)                                     
      370 +    A = zeros(ComplexF64, n, n)                               
      371 +    for (j, q) in enumerate(component)                        
      372 +        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      373 +        for (r, c) in zip(edges.strings, edges.coeffs)        
      374 +            i = get(index, r, 0)                              
      375 +            i == 0 && throw(ArgumentError("orbit component is 
          +not closed"))                                                 
      376 +            A[i, j] += 1im * c / hbar                         
      377 +        end                                                   
      378 +    end                                                       
      379 +    return A                                                  
      380 +end                                                           
      381 +                                                              
      382 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, hbar::Real, componentm
          +ax::Int, krylovdim::Int, truncation)                          
      383 +    componentmax > 0 || throw(ArgumentError("componentmax must
          + be positive"))                                               
      384 +    maxlength = 1000                                          
      385 +    out = zero(O)                                             
      386 +    assigned = Set{eltype(O.strings)}()                       
      387 +    coeff_lookup = Dict{eltype(O.strings),eltype(O.coeffs)}() 
      388 +    for (q, c) in zip(O.strings, O.coeffs)                    
      389 +        coeff_lookup[q] = get(coeff_lookup, q, zero(eltype(O.c
          +oeffs))) + c                                                  
      390 +    end                                                       
      391 +                                                              
      392 +    for q0 in O.strings                                       
      393 +        q0 in assigned && continue                            
      394 +        component, index, capped = _orbit_component(Ha, q0, co
          +mponentmax, maxlength)                                        
      395 +        if capped                                             
      396 +            c0 = get(coeff_lookup, q0, zero(eltype(O.coeffs)))
      397 +            seed_op = typeof(O)([q0], [c0])                   
      398 +            out = out + _orbit_flow_krylov(Ha, seed_op, dt, hb
          +ar, krylovdim, truncation)                                    
      399 +            push!(assigned, q0)                               
      400 +            continue                                          
      401 +        end                                                   
      402 +                                                              
      403 +        A = _orbit_component_matrix(Ha, component, index, hbar
          +, maxlength)                                                  
      404 +        coeffs = zeros(ComplexF64, length(component))         
      405 +        has_active = false                                    
      406 +        for (q, c) in coeff_lookup                            
      407 +            j = get(index, q, 0)                              
      408 +            if j != 0                                         
      409 +                coeffs[j] += c                                
      410 +                push!(assigned, q)                            
      411 +                has_active = true                             
      412 +            end                                               
      413 +        end                                                   
      414 +        has_active || continue                                
      415 +        coeffs2 = exp(dt * A) * coeffs                        
      416 +        out = out + typeof(O)(component, coeffs2)             
      417 +    end                                                       
      418 +    return truncation(out)                                    
      419 +end                                                           
      420 +                                                              
      421  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
           runcation)
      351 -    return _orbit_flow_krylov(Ha, O, dt, hbar, method.krylovdi
          -m, truncation)                                                
      422 +    return _orbit_flow_component(Ha, O, dt, hbar, method.compo
          +nentmax, method.krylovdim, truncation)                        
      423  end
      424  
      425  function _orbit_terms(H::Operator{<:PauliStringTS})

  Ran 2 shell commands (ctrl+o to expand)

● The benchmark is running; I’ll report the results when it completes.

  Read 1 file (ctrl+o to expand)

● Task Output(non-blocking) be7ro42xu
  ⎿  0.0%┣                                                ┫ 0/1 [00:00<00:0
     0, -0s/it]
     100.0%┣██████████████████████████████████████████┫ 1/1 [00:00<Inf:Inf,
     … +251 lines (ctrl+o to expand)

● Stop Task
  ⎿  julia --project=/home/dhruva/Programs/unitaryHack/PauliStrings.jl
     /home/dhruva/.claude/jobs/513a65b3/tmp/compare_active_edges_rk4.jl ·
     stopped

● Implemented the first exact component-exponential backend into TrotterTS.

  What changed

  TrotterTS now defaults to component exponentials with Krylov fallback:

  TrotterTS(; order=2, krylovdim=4, componentmax=1024)

  The orbit flow now does:

  1. Build the connected Pauli-orbit component from orbit_edges.
  2. Build the finite Liouvillian matrix on that component.
  3. Apply the exact dense component exponential:

  exp(dt * A) * coeffs
  4. If the component exceeds componentmax, fall back to Krylov for that seed.

  So this is now genuinely closer to the TS Pauli rotation rule.

  Test status

  test/evolution.jl passes:

  evolve | 8 passed

  Important caveat

  A benchmark on larger cases became extremely slow, so I stopped it. That was
  expected for this first implementation because it still:

  - rebuilds components every orbit flow,
  - rebuilds dense matrices every orbit flow,
  - computes many dense exponentials,
  - has no component plan cache.

  This implementation is correct-oriented, not performant yet.

  Why this is still progress

  We now have the actual orbit-rotation structure in the code:

  ordinary Pauli rotation:
      exact exp on a 2-node Pauli component

  TS Pauli-orbit rotation:
      exact exp on a connected Pauli-orbit component

  The next step is not another algorithm change; it is to cache component plans
  so repeated flows do not rebuild the graph and matrix each step.

✻ Churned for 12m 45s

❯ go ahead                                                                      

● I’ll add per-orbit component plan caches: build each connected component and
  Liouvillian matrix once, cache it by all nodes in the component, and cache
  exp(dt*A) per step size. Then I’ll benchmark small cases to verify the rebuild
  cost actually drops.

  Read 1 file (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 47 lines, removed 10 lines
      347      return truncation(_axpy_combination(V, coeffs, O))
      348  end
      349  
      350 +struct _OrbitComponentPlan{P}                                 
      351 +    component::Vector{P}                                      
      352 +    index::Dict{P,Int}                                        
      353 +    matrix::Matrix{ComplexF64}                                
      354 +    exp_cache::Dict{Float64,Matrix{ComplexF64}}               
      355 +    capped::Bool                                              
      356 +end                                                           
      357 +                                                              
      358 +mutable struct _OrbitFlowCache{P,O}                           
      359 +    plans::Dict{P,_OrbitComponentPlan{P}}                     
      360 +    fallback::Dict{P,O}                                       
      361 +end                                                           
      362 +                                                              
      363 +_OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlow
          +Cache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())     
      364 +                                                              
      365  function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
           :PauliStringTS, componentmax::Int, maxlength::Int)
      366      component = typeof(seed)[seed]
      367      index = Dict{typeof(seed),Int}(seed => 1)
     ...
      394      return A
      395  end
      396  
      382 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, hbar::Real, componentm
          -ax::Int, krylovdim::Int, truncation)                          
      397 +function _component_plan(cache::_OrbitFlowCache, Ha::Operator{
          +<:PauliStringTS}, seed::PauliStringTS, hbar::Real, componentma
          +x::Int, maxlength::Int)                                       
      398 +    plan = get(cache.plans, seed, nothing)                    
      399 +    plan !== nothing && return plan                           
      400 +                                                              
      401 +    component, index, capped = _orbit_component(Ha, seed, comp
          +onentmax, maxlength)                                          
      402 +    A = capped ? zeros(ComplexF64, 0, 0) : _orbit_component_ma
          +trix(Ha, component, index, hbar, maxlength)                   
      403 +    plan = _OrbitComponentPlan(component, index, A, Dict{Float
          +64,Matrix{ComplexF64}}(), capped)                             
      404 +    for q in component                                        
      405 +        cache.plans[q] = plan                                 
      406 +    end                                                       
      407 +    return plan                                               
      408 +end                                                           
      409 +                                                              
      410 +function _component_exp(plan::_OrbitComponentPlan, dt::Real)  
      411 +    key = Float64(dt)                                         
      412 +    E = get(plan.exp_cache, key, nothing)                     
      413 +    if E === nothing                                          
      414 +        E = exp(key * plan.matrix)                            
      415 +        plan.exp_cache[key] = E                               
      416 +    end                                                       
      417 +    return E                                                  
      418 +end                                                           
      419 +                                                              
      420 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          +, hbar::Real, componentmax::Int, krylovdim::Int, truncation)  
      421      componentmax > 0 || throw(ArgumentError("componentmax must
            be positive"))
      422      maxlength = 1000
      423      out = zero(O)
     ...
      429  
      430      for q0 in O.strings
      431          q0 in assigned && continue
      394 -        component, index, capped = _orbit_component(Ha, q0, co
          -mponentmax, maxlength)                                        
      395 -        if capped                                             
      432 +        plan = _component_plan(cache, Ha, q0, hbar, componentm
          +ax, maxlength)                                                
      433 +        if plan.capped                                        
      434              c0 = get(coeff_lookup, q0, zero(eltype(O.coeffs)))
      435              seed_op = typeof(O)([q0], [c0])
      436              out = out + _orbit_flow_krylov(Ha, seed_op, dt, hb
           ar, krylovdim, truncation)
     ...
      438              continue
      439          end
      440  
      403 -        A = _orbit_component_matrix(Ha, component, index, hbar
          -, maxlength)                                                  
      404 -        coeffs = zeros(ComplexF64, length(component))         
      441 +        coeffs = zeros(ComplexF64, length(plan.component))    
      442          has_active = false
      443          for (q, c) in coeff_lookup
      407 -            j = get(index, q, 0)                              
      444 +            j = get(plan.index, q, 0)                         
      445              if j != 0
      446                  coeffs[j] += c
      447                  push!(assigned, q)
     ...
      449              end
      450          end
      451          has_active || continue
      415 -        coeffs2 = exp(dt * A) * coeffs                        
      416 -        out = out + typeof(O)(component, coeffs2)             
      452 +        coeffs2 = _component_exp(plan, dt) * coeffs           
      453 +        out = out + typeof(O)(plan.component, coeffs2)        
      454      end
      455      return truncation(out)
      456  end
      457  
      421 -function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
          -r{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, t
          -runcation)                                                    
      422 -    return _orbit_flow_component(Ha, O, dt, hbar, method.compo
          -nentmax, method.krylovdim, truncation)                        
      458 +function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
          +r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
          +TrotterTS, hbar::Real, truncation)                            
      459 +    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          +d.componentmax, method.krylovdim, truncation)                 
      460  end
      461  
      462  function _orbit_terms(H::Operator{<:PauliStringTS})
● Update(src/evolution.jl)
  ⎿  Added 6 lines, removed 6 lines
      468      return terms
      469  end
      470  
      471 -function _trotterts_step(Hterms, O::Operator{<:PauliStringTS},
          - dt::Real, method::TrotterTS, hbar::Real, truncation)         
      471 +function _trotterts_step(Hterms, caches, O::Operator{<:PauliSt
          +ringTS}, dt::Real, method::TrotterTS, hbar::Real, truncation) 
      472      method.order ∈ (1, 2) || throw(ArgumentError("order must b
           e 1 or 2, got $(method.order)"))
      473      isempty(Hterms) && return O
      474      if method.order == 1 || length(Hterms) == 1
      475 -        for Ha in Hterms                                      
      476 -            O = _orbit_flow(Ha, O, dt, method, hbar, truncatio
          -n)                                                            
      475 +        for j in eachindex(Hterms)                            
      476 +            O = _orbit_flow(Hterms[j], O, dt, caches[j], metho
          +d, hbar, truncation)                                          
      477          end
      478      else
      479          L = length(Hterms)
      480          for j in 1:(L - 1)
      481 -            O = _orbit_flow(Hterms[j], O, dt / 2, method, hbar
          -, truncation)                                                 
      481 +            O = _orbit_flow(Hterms[j], O, dt / 2, caches[j], m
          +ethod, hbar, truncation)                                      
      482          end
      483 -        O = _orbit_flow(Hterms[L], O, dt, method, hbar, trunca
          -tion)                                                         
      483 +        O = _orbit_flow(Hterms[L], O, dt, caches[L], method, h
          +bar, truncation)                                              
      484          for j in (L - 1):-1:1
      485 -            O = _orbit_flow(Hterms[j], O, dt / 2, method, hbar
          -, truncation)                                                 
      485 +            O = _orbit_flow(Hterms[j], O, dt / 2, caches[j], m
          +ethod, hbar, truncation)                                      
      486          end
      487      end
      488      return O

● Update(src/evolution.jl)
  ⎿  Added 2 lines, removed 1 line
      495      n = length(tspan)
      496      history = _alloc_history(fout, O, n)
      497      Hterms = _orbit_terms(H)
      498 +    caches = [_OrbitFlowCache(paulistringtype(H), typeof(H)) f
          +or _ in Hterms]                                               
      499      O = copy(O)
      500      for i in ProgressBar(1:(n - 1))
      501          dt = tspan[i + 1] - tspan[i]
      501 -        O = _trotterts_step(Hterms, O, dt, method, hbar, trunc
          -ation)                                                        
      502 +        O = _trotterts_step(Hterms, caches, O, dt, method, hba
          +r, truncation)                                                
      503          O = dissipation(O, dt)
      504          O = truncation(O)
      505          _save!(history, fout, O, i + 1)

  Ran 3 shell commands (ctrl+o to expand)

● Implemented component plan caching, but performance is still bad.

  What changed

  TrotterTS now builds per-orbit caches:

  _OrbitFlowCache

  Each cache stores component plans:

  _OrbitComponentPlan(
      component,
      index,
      matrix,
      exp_cache,
      capped
  )

  When a component is discovered, the plan is registered for every node in that
  component. Matrix exponentials are cached by dt, so repeated half/full steps
  reuse:

  exp(dt * A)

  This avoids rebuilding the same component/matrix/exponential repeatedly within
  one evolution.

  Tests

  test/evolution.jl passes.

  Benchmark

  N=16, M=10, steps=20

  ┌────────────────────────────┬─────────┬──────────┬───────────┐
  │           Method           │  Time   │  Alloc   │ Final len │
  ├────────────────────────────┼─────────┼──────────┼───────────┤
  │ RK4                        │ 0.0569s │ 0.073 GB │      1024 │
  ├────────────────────────────┼─────────┼──────────┼───────────┤
  │ Expanded Trotter           │ 0.0697s │ 0.116 GB │       514 │
  ├────────────────────────────┼─────────┼──────────┼───────────┤
  │ TrotterTS component cached │   2.11s │  3.67 GB │      1024 │
  └────────────────────────────┴─────────┴──────────┴───────────┘

  N=32, M=12, steps=10

  ┌────────────────────────────┬────────┬──────────┬───────────┐
  │           Method           │  Time  │  Alloc   │ Final len │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ RK4                        │ 0.117s │ 0.102 GB │      4096 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ Expanded Trotter           │ 0.236s │ 0.412 GB │      1401 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ TrotterTS component cached │ 11.07s │ 20.09 GB │      4096 │
  └────────────────────────────┴────────┴──────────┴───────────┘

  Interpretation

  Caching components/matrices/exponentials is not enough. The component approach
  is still far too allocation-heavy.

  Likely causes:

  1. Many components are built and cached as full Vector + Dict + dense Matrix.
  2. Even small dense exponentials/matrix-vector multiplies repeated across
  hundreds of components create substantial overhead.
  3. Each orbit flow still gathers/scatters by scanning coefficient dictionaries
  and component dictionaries.
  4. The active basis changes under truncation, so caching helps but does not
  eliminate component discovery.

  Conclusion

  The exact component-exponential rule is mathematically correct, but this
  direct implementation is not viable.

  At this point, the core theory is sound, but the practical backend remains
  unresolved:

  - Krylov is conceptually clean but too much repeated high-level OperatorTS
  work.
  - Component exponentials are exact but too much graph/matrix bookkeeping.
  - orbit_edges shows the graph is sparse, but using it naively creates
  overhead.

  I would not continue optimizing this component backend without redesigning the
  data representation around compact indexed orbit graph components.

✻ Cogitated for 3m 41s

❯ ok go ahead                                                                   

● I’ll first fix the obvious component-backend inefficiency: it currently
  gathers coefficients by scanning the entire active operator for every
  component. That is O(#components × length(O)) and explains the terrible
  runtime. I’ll change it to gather by iterating each component’s own basis
  only, then benchmark again.

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 3 lines
      440  
      441          coeffs = zeros(ComplexF64, length(plan.component))
      442          has_active = false
      443 -        for (q, c) in coeff_lookup                            
      444 -            j = get(plan.index, q, 0)                         
      445 -            if j != 0                                         
      443 +        for (j, q) in enumerate(plan.component)               
      444 +            c = get(coeff_lookup, q, nothing)                 
      445 +            if c !== nothing                                  
      446                  coeffs[j] += c
      447                  push!(assigned, q)
      448                  has_active = true

  Ran 1 shell command (ctrl+o to expand)

● Added component plan caches and fixed the obvious O(#components × length(O))
  gather issue.

  What changed

  For each orbit term, TrotterTS now maintains an _OrbitFlowCache containing:

  Dict{PauliStringTS, _OrbitComponentPlan}

  Each _OrbitComponentPlan stores:

  - component basis,
  - index map,
  - Liouvillian matrix,
  - cached matrix exponentials keyed by dt,
  - capped/fallback flag.

  When a component is built, the same plan is registered for every node in the
  component, so future active strings reuse it.

  Also optimized coefficient gather:

  Before, for each component it scanned the entire active operator:

  for each component:
      for each q in O:
          check if q in component

  Now it iterates only over the component basis:

  for each q in component:
      check coefficient lookup

  Tests

  test/evolution.jl passes.

  Benchmark after caching and gather fix

  N=16, M=10, steps=20

  ┌────────────────────────────┬────────┬──────────┬───────────┐
  │           Method           │  Time  │  Alloc   │ Final len │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ RK4                        │ 0.058s │ 0.076 GB │      1024 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ Expanded Trotter           │ 0.073s │ 0.116 GB │       514 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ TrotterTS cached component │ 1.425s │ 3.718 GB │      1024 │
  └────────────────────────────┴────────┴──────────┴───────────┘

  Before gather fix, same rough case was:

  TrotterTS component: ~2.11s

  So gather fix helped, but not enough.

  N=32, M=12, steps=5

  ┌────────────────────────────┬────────┬──────────┬───────────┐
  │           Method           │  Time  │  Alloc   │ Final len │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ RK4                        │ 0.026s │ 0.024 GB │      4096 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ Expanded Trotter           │ 0.090s │ 0.169 GB │      1221 │
  ├────────────────────────────┼────────┼──────────┼───────────┤
  │ TrotterTS cached component │ 1.978s │ 5.582 GB │      4096 │
  └────────────────────────────┴────────┴──────────┴───────────┘

  Interpretation

  Caching plans and fixing coefficient gather reduced some waste, but component
  exponentials are still far too expensive.

  The dominant issue is likely:

  many components × dense matrix construction/exp/matvec/scatter overhead

  Even if components are small, there are hundreds of them.

  At N=32, earlier stats showed ~700 components per orbit. With 5 orbit factors
  per step, that is thousands of component operations.

  So the component-exponential backend is exact and conceptually right, but the
  naive “one dense exp per component” strategy is not viable without more
  aggressive global planning or batching.

  Current state of candidate backends

  - Krylov: conceptually correct, too many OperatorTS operations.
  - Component exact exp: mathematically ideal, too many tiny component
  plans/exps.
  - Taylor2: fast-ish but not the right Trotter primitive.
  - Expanded Trotter: exact product over ordinary translated gates, wrong
  granularity.

  This suggests the next conceptual step is probably not per-component dense
  exponentials, but a way to batch many identical/isomorphic components or
  derive closed forms for common component types — while staying generic if
  possible.

✻ Crunched for 1m 27s

❯ what are these components?                                                    

● In this context, a component is a connected component of the Pauli-orbit graph
  generated by one Hamiltonian orbit.

  Let one TS Hamiltonian orbit be:

  H_a^{orb} = c_a \sum_r T_r P_a T_r^\dagger.

  This defines a Liouvillian:

  L_a(O) = \frac{i}{\hbar}[H_a^{orb}, O].

  Now take one TS Pauli orbit [Q]. Applying L_a to it gives a linear combination
  of other TS Pauli orbits:

  L_a([Q]) = \sum_R A_{R,Q} [R].

  Those nonzero transitions are graph edges:

  [Q] -> [R]

  where [R] comes from commuting Q with translated copies of P_a, then folding
  back into canonical PauliStringTS representatives.

  A component is the set of all TS Pauli orbits reachable from [Q] by repeatedly
  applying those edges.

  So:

  component([Q]) =
  { [Q],
    all orbits reachable from L_a[Q],
    all orbits reachable from L_a²[Q],
    all orbits reachable from L_a³[Q],
    ... }

  Example analogy

  For ordinary non-TS Pauli rotation with generator P:

  e^{i\theta P/2} Q e^{-i\theta P/2}

  the graph component is tiny:

  Q <-> P Q

  So the component is usually size 1 or 2.

  That is why ordinary Pauli rotation has the closed form:

  Q \mapsto \cos(\theta) Q + \sin(\theta) Q'.

  For a TS orbit generator [P], the component can be larger:

  [Q] -> [R1], [R2], [R3], ...

  but it may still be small enough to exponentiate exactly.

  Concrete meaning in code

  If:

  q::PauliStringTS
  Ha::OperatorTS   # one orbit term

  then:

  orbit_edges(Ha, q)

  returns all immediate neighbors of q:

  q -> r₁, r₂, ...

  The component builder repeatedly calls orbit_edges:

  component = [q]
  frontier = [q]

  while frontier not empty
      for q in frontier
          edges = orbit_edges(Ha, q)
          add unseen edge targets to component
      end
  end

  So the component is a finite basis of TS Pauli strings on which L_a closes.

  Why components matter

  On one component, the orbit Liouvillian becomes a finite matrix:

  A_C

  Then the exact orbit rotation is:

  e^{\tau A_C}.

  So the component is the TS analogue of the two-dimensional {Q, PQ} subspace in
  an ordinary Pauli rotation.

  What we observed

  For sampled evolved states, many components were small:

  p50 ~ 4–8
  p90 ~ 16–32

  But there were many components. That is why the naive implementation was slow:

  hundreds of small components × dense exp per component

  Even though each component is small, the overhead of building/exponentiating
  many of them dominates.

✻ Cooked for 26s

❯ I don;t know what to do next then. If number of components are high, then any 
  speedup in the exp per component will not translate to genuine speedups.      

● I think that concern is correct.

  If an orbit flow decomposes into hundreds or thousands of small components,
  then even if each component exponential is cheap, the overhead of:

  component lookup
  coefficient gather
  small matrix-vector multiply
  scatter
  bookkeeping

  can dominate. This is exactly what we saw.

  So component exponentials are mathematically clean, but probably not the
  practical route unless we can batch many identical components. That may be
  possible, but it is a much more specialized representation problem.

  Where we stand

  We have now learned:

  1. Expanded TS Trotter is slow because it applies O(Nm) ordinary gates.
  2. Orbit-level Trotter is the right mathematical grouping.
  3. Taylor-2 orbit flow is too RK-like, so not the desired solution.
  4. Krylov orbit flow is too expensive because it repeatedly builds Krylov
  vectors as OperatorTS.
  5. Exact component exponentials are too granular because there are many
  components.
  6. Generic active orbit edges show graph sparsity, but using the graph
  directly has overhead.

  So the missing piece is not just implementation polish. It is choosing the
  right finite orbit-flow representation.

  Possible next directions

  Option A: Batch isomorphic components

  Many small components may have the same graph structure and differ only by
  translated/canonical labels.

  If true, then we can group components by a graph signature:

  same component matrix A
  different coefficient vectors

  Then one cached exponential applies to many components.

  This could convert:

  hundreds of exp(A_i)

  into:

  a few exp(A_signature)
  many small matvecs

  But this requires proving or measuring that component types repeat heavily.

  Next check: compute component signatures and count unique signatures.

  If there are few signatures, component method is back on the table.
  If almost every component is unique, abandon this route.

  ---
  Option B: Orbit-flow as direct sparse matrix on active basis

  Instead of components, build one sparse matrix for L_a on the active/truncated
  basis:

  A_a

  Then apply:

  e^{τ A_a}

  approximately using sparse exponential action or maybe low-degree
  rational/Chebyshev.

  But this starts looking like Krylov again, only with a better data structure.

  It could be worthwhile if we build A_a once and reuse it many steps, but
  truncation changes the active basis.

  ---
  Option C: Accept TS-native Trotter needs extra structure

  This may be the uncomfortable theoretical conclusion:

  ▎ Translation symmetry alone gives representation compression, but not 
  ▎ necessarily a cheap finite-time orbit exponential.

  The non-TS Trotter advantage comes from exact single-Pauli rotations. Once we
  group translated Paulis into an orbit, generic exact orbit flows may be too
  complex unless the orbit has extra structure.

  Extra structures include:

  - commuting translated copies,
  - graph colorability,
  - low component type count,
  - locality with reusable patterns,
  - fixed active basis.

  But “translation symmetry alone” may not guarantee a fast Trotter primitive.

  This does not mean failure; it means the paper/PR should be honest: TS-native
  orbit Trotter is mathematically right, but generic speedup requires an
  efficient orbit-flow backend, and that may need additional exploitable
  structure.

  ---
  Option D: Go back to current expanded Trotter but improve truncation semantics

  This is pragmatic, not theoretically pure.

  Use expanded gates but keep TS compression and tune truncation/folding. We
  already found this is only a partial fix. It will not be the big algorithmic
  win, but it might improve package behavior.

  ---
  My recommendation

  Before giving up on component exponentials, do one decisive check:

  ▎ Count unique component signatures.

  For each orbit term and a representative evolved state:

  1. Build components.
  2. For each component, compute a canonical signature of its internal graph:
    - component size,
    - sorted weighted adjacency pattern up to relabeling? initially maybe
  simpler:
        - size,
      - multiset of row/column degrees,
      - sorted eigenvalues of component matrix rounded.
  3. Count:
    - number of components,
    - number of unique signatures,
    - frequency distribution.

  If we see:

  700 components but 5 signatures

  then batching can work.

  If we see:

  700 components and 650 signatures

  then component exponentials are dead as a generic route.

  This is a small, informative experiment and directly addresses your concern.

✻ Worked for 29s

❯ Ok do it. Let me again caution: if we are generating components from the TS   
  state, then aren't we losing the TS representation? I understand that the     
  graph is sparse, but still. anyway, continue                                  

● That caution is valid. The components are still built in the TS Pauli-orbit 
  basis (PauliStringTS nodes), not in ordinary resum space, so we are not
  literally losing TS representation. But we may lose the compact active 
  representation if each active TS string belongs to a distinct component with
  its own bookkeeping. That’s exactly why this signature check matters: many
  components can erase the practical benefit even though every node is
  TS-compressed.

● Write(~/.claude/jobs/513a65b3/tmp/component_signature_stats.jl)
  ⎿  Wrote 122 lines to
     ../../../.claude/jobs/513a65b3/tmp/component_signature_stats.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
       8 function Xtot(N)
       9     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
      10 end
     … +123 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I measured component signature multiplicities.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/component_signature_stats.jl

  This builds components in the TS Pauli-orbit basis, computes:

  - component size distribution,
  - spectral signature of the component matrix,
  - cheaper degree signature,
  - multiplicities of repeated signatures.

  Important clarification

  We are not expanding into ordinary Pauli strings. Components are sets of
  PauliStringTS nodes.

  But your concern is still valid: if the TS state decomposes into
  hundreds/thousands of separate components, the component bookkeeping can
  destroy the practical benefit of TS compression.

  Results

  N=16, M=10, lenO=1024

  ┌───────┬────────────┬─────────────────┬───────────────┬─────┬─────┬─────┐
  │ Orbit │ Components │ Unique spectral │ Unique degree │ p50 │ p90 │ Max │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     1 │        251 │              40 │             7 │   8 │  16 │  64 │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     2 │        217 │              25 │             4 │   4 │  16 │  64 │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     3 │        229 │              50 │             8 │   4 │  32 │ 128 │
  └───────┴────────────┴─────────────────┴───────────────┴─────┴─────┴─────┘

  Top spectral multiplicities:

  orbit 1: [69, 39, 37, 36, 11, ...]
  orbit 2: [150, 11, 8, 5, 3, ...]
  orbit 3: [64, 40, 30, 27, 11, ...]

  N=32, M=12, lenO=4096

  ┌───────┬────────────┬─────────────────┬───────────────┬─────┬─────┬─────┐
  │ Orbit │ Components │ Unique spectral │ Unique degree │ p50 │ p90 │ Max │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     1 │        782 │             121 │             8 │   8 │  32 │ 128 │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     2 │        670 │             129 │             9 │   8 │  32 │ 256 │
  ├───────┼────────────┼─────────────────┼───────────────┼─────┼─────┼─────┤
  │     3 │        697 │              44 │             4 │   4 │  16 │  64 │
  └───────┴────────────┴─────────────────┴───────────────┴─────┴─────┴─────┘

  Top spectral multiplicities:

  orbit 1: [179, 112, 109, 92, 20, ...]
  orbit 2: [152, 90, 86, 71, 19, ...]
  orbit 3: [398, 19, 19, 16, 14, ...]

  N=64, M=14, sampled first 300 components per orbit

  Partial active assignment only because of sample cap.

  ┌───────┬─────────────────┬───────────────┬─────────────┬─────┬─────┬─────┐
  │ Orbit │    Sampled      │    Unique     │   Unique    │ p50 │ p90 │ Max │
  │       │   components    │   spectral    │   degree    │     │     │     │
  ├───────┼─────────────────┼───────────────┼─────────────┼─────┼─────┼─────┤
  │     1 │             300 │            58 │           8 │   8 │  32 │ 128 │
  ├───────┼─────────────────┼───────────────┼─────────────┼─────┼─────┼─────┤
  │     2 │             300 │            49 │           7 │   8 │  16 │  64 │
  ├───────┼─────────────────┼───────────────┼─────────────┼─────┼─────┼─────┤
  │     3 │             300 │            27 │           4 │   4 │  16 │  64 │
  └───────┴─────────────────┴───────────────┴─────────────┴─────┴─────┴─────┘

  Top spectral multiplicities:

  orbit 1: [77, 50, 42, 38, 14, ...]
  orbit 2: [84, 49, 43, 37, 11, ...]
  orbit 3: [197, 14, 8, 7, 6, ...]

  Interpretation

  There is repetition, but not enough to trivially fix the component backend.

  Good news:

  - degree signatures are few: ~4–9 unique in these cases.
  - some spectral signatures have high multiplicity, especially orbit 3.
  - components are generally small.

  Bad news:

  - spectral signatures still number dozens to >100 at N=32.
  - component count is high: ~700 components per orbit at N=32.
  - even if we reuse exponentials by signature, we still need gather/scatter for
  many components.

  What this means

  Batching identical components may help, but it is not an obvious silver
  bullet.

  The best case is orbit 3 at N=32:

  697 components, 44 spectral signatures, top one appears 398 times

  That is promising.

  But orbit 1 and 2 are less clean:

  ~700 components, ~120 spectral signatures

  Still maybe batchable, but more involved.

  Conclusion

  Your concern is confirmed:

  ▎ Many components can erase the practical benefit, even though components live
  ▎ in TS representation.

  Component exponentials remain theoretically valid, but practical speedup
  likely requires more than caching individual components. We would need
  batching by component type/signature or a more global indexed representation.

  At this point, the most defensible path is to document that the exact TS
  orbit-rotation rule exists, but that naive component decomposition is not yet
  a performant backend.

✻ Worked for 1m 23s

❯ Is it possible to drop some components, on the basis of their weights? So the 
  trimming is happening as part of the update, not as a postprocess step        

● Yes, that is possible, and it may be exactly the right way to make component
  propagation practical.

  Instead of:

  propagate all components exactly
  then trim the final OperatorTS

  we can do:

  estimate component importance
  skip or approximate low-weight components
  propagate important components exactly

  This makes truncation part of the orbit update.

  Why this is natural

  The current truncation policy trim(O, M) keeps the largest individual
  Pauli-orbit coefficients after an update.

  But component propagation gives a coarser unit:

  component C = set of TS Pauli orbits coupled by one orbit Liouvillian

  During an orbit flow, components do not mix. So each component carries its own
  coefficient vector:

  c_C.

  The orbit flow is:

  c_C' = e^{\tau A_C} c_C.

  A natural component weight is:

  w_C = \|c_C\|^2

  or Hilbert-Schmidt weighted norm:

  w_C = \sum_{q \in C} |c_q|^2 \|[q]\|^2.

  Then we can keep only the most important components before applying expensive
  exact exponentials.

  Component-aware truncation rule

  For one orbit flow:

  1. Decompose active O into components under L_a.
  2. For each component C, compute active coefficient weight:

  w_C = \sum_{q \in C \cap supp(O)} |c_q|^2

  2. or metric-weighted version.
  3. Select components until budget is reached:
    - keep top K components by weight,
    - or keep components whose cumulative weight exceeds threshold,
    - or keep components covering expected M output strings.
  4. For kept components:
    - apply exact component exponential.
  5. For dropped components:
    - either discard them,
    - or leave them unchanged,
    - or use cheaper first-order/Taylor/Krylov-low-dim approximation.
  6. Then optionally run normal trim after the orbit flow as a safety step.

  Drop versus freeze

  There are two choices.

  A. Drop low-weight components

  c_C' = 0

  Pros:
  - aggressive truncation,
  - fewer terms,
  - fastest.

  Cons:
  - may lose important conserved/norm contributions,
  - can be too dissipative.

  B. Freeze low-weight components

  c_C' = c_C

  Pros:
  - does not delete information immediately,
  - cheap,
  - more stable.

  Cons:
  - inaccurate for those components,
  - may keep too many terms.

  C. Cheap approximate update

  For low-weight components:

  c_C' \approx c_C + \tau A_C c_C

  or no component matrix; just one orbit_edges application.

  Pros:
  - middle ground.

  Cons:
  - reintroduces RK-like updates.

  I would test both drop and freeze. Dropping is closer to integrated trimming;
  freezing is safer.

  Why this could help

  Previously component exact propagation was slow because we had:

  hundreds of components per orbit

  But coefficient weight is usually heavy-tailed. If most norm is in a small
  fraction of components, then we only exponentiate those.

  This changes cost from:

  all components

  to:

  important components only

  This is analogous to trim(O, M), but at the component level.

  This also addresses your TS concern

  We remain in OperatorTS, but instead of tracking every TS component, we choose
  which orbit components matter. That is a TS-native truncation strategy.

  So the update becomes:

  TS orbit flow + TS component-aware truncation

  rather than:

  full exact orbit flow + post hoc coefficient trim

  What to test

  Before implementing into core, write a prototype that reports component
  weights.

  For an evolved O, for each orbit Ha:

  number of components
  component sizes
  component active weights
  cumulative weight distribution
  components needed for 90%, 99%, 99.9% of weight
  number of active strings covered

  If, for example:

  top 10% components carry 99% of weight

  then this is promising.

  If weights are flat across hundreds of components, it won’t help.

  Proposed prototype

  For each orbit flow:

  components = build_components(Ha, O)
  weights = component_weight.(components, O)

  sort by weight descending

  keep components while:
      cumulative_weight / total_weight < threshold
      or output_size_budget not exceeded

  Then compare:

  - exact full component flow,
  - component-truncated flow,
  - RK4,
  - expanded Trotter,

  on observable accuracy and runtime.

  Important: metric-weighted component norm

  Because OperatorTS orbit representatives may have different orbit
  multiplicities/stabilizers, the best weight is probably:

  trace_product(O_C', O_C)

  or equivalent metric-weighted norm, not just coefficient norm.

  But for a first diagnostic, plain coefficient norm is fine.

  Bottom line

  Yes: we can make trimming part of the TS orbit update by dropping/freezing
  low-weight components before exact component exponentiation.

  This is promising because it attacks the main component-backend failure mode:

  too many components

  without abandoning the TS-native orbit-rotation rule.

✻ Churned for 37s

❯ /model                                                                        
  ⎿  Set model to openai.gpt-5.5 and saved as your default for new sessions

❯ So find a method to reduce components then, possibly revisit the math         

● I’ll first test whether component weights are heavy-tailed. If most weight
  sits in a small fraction of components, then component-level truncation is
  mathematically plausible. If weights are flat, reducing components will likely
  damage accuracy.

● Write(~/.claude/jobs/513a65b3/tmp/component_weight_distribution.jl)
  ⎿  Wrote 127 lines to
     ../../../.claude/jobs/513a65b3/tmp/component_weight_distribution.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function grow(H,O0,M,steps)
     … +128 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)
● Update(~/.claude/jobs/513a65b3/tmp/component_weight_distribution.jl)
  ⎿  Added 1 line, removed 1 line
      83          w,a = use_metric ? metric_weight(comp, coeffs, typeof(O
          )) : plain_weight(comp, coeffs)
      84          push!(weights,w); push!(sizes,length(comp)); push!(acti
          ves,a)
      85          for q in comp
      86 -            q in coeffs && push!(assigned,q)                   
      86 +            haskey(coeffs, q) && push!(assigned,q)             
      87          end
      88      end
      89      order=sortperm(weights; rev=true)

  Ran 1 shell command (ctrl+o to expand)

● I measured component weight concentration to find a way to reduce components
  as part of the orbit update.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/component_weight_distribution.jl

  Key result

  Component weights are extremely heavy-tailed.

  For each orbit flow, only a tiny fraction of components carry almost all the
  coefficient weight.

  Plain coefficient weights

  N=16, M=10, lenO=1024

  ┌───────┬──────────────┬────────────┬─────┬───────┬────────┐
  │ Orbit │ # components │ 90% weight │ 99% │ 99.9% │ 99.99% │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     1 │          251 │    4 comps │   8 │    19 │     36 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     2 │          217 │    5 comps │  12 │    19 │     30 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     3 │          229 │    5 comps │  13 │    24 │     38 │
  └───────┴──────────────┴────────────┴─────┴───────┴────────┘

  N=32, M=12, lenO=4096

  ┌───────┬──────────────┬────────────┬─────┬───────┬────────┐
  │ Orbit │ # components │ 90% weight │ 99% │ 99.9% │ 99.99% │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     1 │          782 │    4 comps │   8 │    19 │     36 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     2 │          670 │    5 comps │  13 │    24 │     38 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     3 │          697 │    5 comps │  12 │    19 │     30 │
  └───────┴──────────────┴────────────┴─────┴───────┴────────┘

  N=64, M=14, lenO=16384

  ┌───────┬──────────────┬────────────┬─────┬───────┬────────┐
  │ Orbit │ # components │ 90% weight │ 99% │ 99.9% │ 99.99% │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     1 │         1869 │    5 comps │  13 │    24 │     38 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     2 │         2466 │    4 comps │   8 │    19 │     36 │
  ├───────┼──────────────┼────────────┼─────┼───────┼────────┤
  │     3 │         2168 │    5 comps │  12 │    19 │     30 │
  └───────┴──────────────┴────────────┴─────┴───────┴────────┘

  So even for thousands of components:

  ~30–40 components carry 99.99% of the weight

  Metric-weighted component weights

  I also checked Hilbert-Schmidt-like metric weights for N=16,32.

  They gave essentially the same component rankings and cutoff counts.

  Example N=32:

  ┌───────┬──────────────┬─────┬─────┬───────┬────────┐
  │ Orbit │ # components │ 90% │ 99% │ 99.9% │ 99.99% │
  ├───────┼──────────────┼─────┼─────┼───────┼────────┤
  │     1 │          782 │   4 │   8 │    19 │     36 │
  ├───────┼──────────────┼─────┼─────┼───────┼────────┤
  │     2 │          670 │   5 │  13 │    24 │     38 │
  ├───────┼──────────────┼─────┼─────┼───────┼────────┤
  │     3 │          697 │   5 │  12 │    19 │     30 │
  └───────┴──────────────┴─────┴─────┴───────┴────────┘

  So plain coefficient weight is likely acceptable for first implementation.

  Interpretation

  This is promising.

  The component-exponential backend failed because it propagated every
  component:

  hundreds/thousands of components

  But the actual state weight lives in very few components.

  So we can make truncation part of the orbit update:

  for each orbit flow:
      build/identify active components
      compute component weights
      propagate only top-weight components exactly
      drop or freeze the rest

  This directly addresses the failure mode.

  Proposed component-reduction rule

  For each orbit flow, compute component weights:

  w_C = \sum_{q \in C \cap supp(O)} |c_q|^2

  Then keep components until cumulative weight reaches a threshold:

  component_weight_threshold = 0.9999

  or until a component count cap is reached.

  Options for discarded components:

  1. Drop: remove their coefficients.
  2. Freeze: carry them forward unchanged.

  I recommend testing freeze first because it is less destructive:

  important components: exact orbit exponential
  low-weight components: unchanged

  This is still a TS-native orbit update with integrated truncation.

  Why this may work

  For N=64, keeping 99.99% weight means only about:

  30–40 components per orbit

  instead of:

  ~2000 components per orbit

  That is a massive reduction.

  This is the first result that makes component exponentials look practically
  viable.

✻ Brewed for 2m 1s

❯ go ahead with the freeze method                                               

● I’ll add component-weight truncation into the component orbit flow: keep top
  components until a cumulative weight threshold and freeze the rest unchanged.
  Then I’ll run exact tests and benchmark against RK4.

● Update(src/evolution.jl)
  ⎿  Added 8 lines, removed 5 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, krylovdim=4, componentmax=1024)       
      29 +    TrotterTS(; order=2, krylovdim=4, componentmax=1024, compon
         +enttol=0.9999)                                                 
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34  Liouvillians. Each orbit flow is applied by exponentiating conn
          ected components of
      35 -the Pauli-orbit graph. Components larger than `componentmax` fa
         -ll back to a                                                   
      36 -matrix-free Arnoldi/Krylov exponential action with dimension `k
         -rylovdim`.                                                     
      35 +the Pauli-orbit graph. Components are processed in descending c
         +oefficient weight                                              
      36 +until `componenttol` cumulative weight is reached; lower-weight
         + components are                                                
      37 +frozen for that orbit flow. Components larger than `componentma
         +x` fall back to a                                              
      38 +matrix-free Arnoldi/Krylov exponential action with dimension `k
         +rylovdim` when kept.                                           
      39  """
      40  struct TrotterTS <: AbstractEvolutionMethod
      41      order::Int
      42      krylovdim::Int
      43      componentmax::Int
      44 +    componenttol::Float64                                      
      45  end
      43 -TrotterTS(; order::Integer=2, krylovdim::Integer=4, componentma
         -x::Integer=1024) =                                             
      44 -    TrotterTS(Int(order), Int(krylovdim), Int(componentmax))   
      46 +TrotterTS(; order::Integer=2, krylovdim::Integer=4, componentma
         +x::Integer=1024, componenttol::Real=0.9999) =                  
      47 +    TrotterTS(Int(order), Int(krylovdim), Int(componentmax), Fl
         +oat64(componenttol))                                           
      48  
      49  """
      50      RK4()

● Update(src/evolution.jl)
  ⎿  Added 42 lines, removed 14 lines
      420      return E
      421  end
      422  
      423 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          -, hbar::Real, componentmax::Int, krylovdim::Int, truncation)  
      423 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          +, hbar::Real, componentmax::Int, krylovdim::Int, componenttol:
          +:Float64, truncation)                                         
      424      componentmax > 0 || throw(ArgumentError("componentmax must
            be positive"))
      425 +    0 < componenttol <= 1 || throw(ArgumentError("componenttol
          + must be in (0, 1]"))                                         
      426      maxlength = 1000
      426 -    out = zero(O)                                             
      427 -    assigned = Set{eltype(O.strings)}()                       
      427      coeff_lookup = Dict{eltype(O.strings),eltype(O.coeffs)}()
      428      for (q, c) in zip(O.strings, O.coeffs)
      429          coeff_lookup[q] = get(coeff_lookup, q, zero(eltype(O.c
           oeffs))) + c
      430      end
      431  
      432 +    assigned = Set{eltype(O.strings)}()                       
      433 +    component_data = Tuple{Float64,Any,Vector{ComplexF64},Bool
          +}[]                                                           
      434 +    total_weight = 0.0                                        
      435 +                                                              
      436      for q0 in O.strings
      437          q0 in assigned && continue
      438          plan = _component_plan(cache, Ha, q0, hbar, componentm
           ax, maxlength)
      436 -        if plan.capped                                        
      437 -            c0 = get(coeff_lookup, q0, zero(eltype(O.coeffs)))
      438 -            seed_op = typeof(O)([q0], [c0])                   
      439 -            out = out + _orbit_flow_krylov(Ha, seed_op, dt, hb
          -ar, krylovdim, truncation)                                    
      440 -            push!(assigned, q0)                               
      441 -            continue                                          
      442 -        end                                                   
      443 -                                                              
      439          coeffs = zeros(ComplexF64, length(plan.component))
      440 +        weight = 0.0                                          
      441          has_active = false
      442          for (j, q) in enumerate(plan.component)
      443              c = get(coeff_lookup, q, nothing)
      444              if c !== nothing
      445                  coeffs[j] += c
      446 +                weight += abs2(c)                             
      447                  push!(assigned, q)
      448                  has_active = true
      449              end
      450          end
      451          has_active || continue
      455 -        coeffs2 = _component_exp(plan, dt) * coeffs           
      456 -        out = out + typeof(O)(plan.component, coeffs2)        
      452 +        total_weight += weight                                
      453 +        push!(component_data, (weight, plan, coeffs, plan.capp
          +ed))                                                          
      454      end
      455 +                                                              
      456 +    sort!(component_data; by=x -> x[1], rev=true)             
      457 +    keep_weight = componenttol * total_weight                 
      458 +    accumulated = 0.0                                         
      459 +    out = zero(O)                                             
      460 +                                                              
      461 +    for (weight, plan, coeffs, capped) in component_data      
      462 +        if accumulated < keep_weight                          
      463 +            if capped                                         
      464 +                # Fallback only for kept components that excee
          +d the exact-component cap.                                    
      465 +                seed_strings = eltype(O.strings)[]            
      466 +                seed_coeffs = eltype(O.coeffs)[]              
      467 +                for (j, q) in enumerate(plan.component)       
      468 +                    !iszero(coeffs[j]) || continue            
      469 +                    push!(seed_strings, q)                    
      470 +                    push!(seed_coeffs, coeffs[j])             
      471 +                end                                           
      472 +                out = out + _orbit_flow_krylov(Ha, typeof(O)(s
          +eed_strings, seed_coeffs), dt, hbar, krylovdim, truncation)   
      473 +            else                                              
      474 +                coeffs2 = _component_exp(plan, dt) * coeffs   
      475 +                out = out + typeof(O)(plan.component, coeffs2)
      476 +            end                                               
      477 +            accumulated += weight                             
      478 +        else                                                  
      479 +            # Freeze low-weight components: carry their active
          + coefficients forward unchanged.                              
      480 +            for (j, q) in enumerate(plan.component)           
      481 +                !iszero(coeffs[j]) || continue                
      482 +                out = out + typeof(O)([q], [coeffs[j]])       
      483 +            end                                               
      484 +        end                                                   
      485 +    end                                                       
      486      return truncation(out)
      487  end
      488  
      489  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
           TrotterTS, hbar::Real, truncation)
      462 -    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          -d.componentmax, method.krylovdim, truncation)                 
      490 +    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          +d.componentmax, method.krylovdim, method.componenttol, truncat
          +ion)                                                          
      491  end
      492  
      493  function _orbit_terms(H::Operator{<:PauliStringTS})

  Ran 1 shell command (ctrl+o to expand)

● Implemented component-weight truncation with the freeze method.

  What changed

  TrotterTS now has:

  componenttol::Float64

  Constructor:

  TrotterTS(; order=2, krylovdim=4, componentmax=1024, componenttol=0.9999)

  During each orbit flow:

  1. Build active components.
  2. Compute component weights:

  weight = sum(abs2(coeffs in component))
  3. Sort components by descending weight.
  4. Propagate components exactly until cumulative weight reaches:

  componenttol * total_weight
  5. Freeze the rest:

  low-weight component coefficients are carried forward unchanged

  If a kept component exceeds componentmax, it falls back to Krylov.

  Test status

  test/evolution.jl passes:

  evolve | 8 passed

  The test took longer, but still passed.

  Current default

  Default is conservative:

  componenttol = 0.9999

  Based on the component-weight measurements, this should often keep only ~30–40
  components per orbit while retaining 99.99% component weight.

  Important caveat

  This is an approximation integrated into the orbit update. It is no longer the
  exact orbit component exponential for all components; low-weight components
  are frozen. This is intentional and corresponds to component-level truncation
  inside the orbit flow.

✻ Crunched for 2m 0s

❯ Ok then profile it                                                            

● Write(~/.claude/jobs/513a65b3/tmp/profile_component_freeze.jl)
  ⎿  Wrote 29 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_component_freeze.jl
       1 using PauliStrings
       2 using Printf
       3 
       4 function MFIM(N,h)
       5     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
       6 end
       7 function Xtot(N)
       8     H=Operator(N); H += "X",1; OperatorTS{(N,)}(H)
       9 end
      10 function bench(N,M,steps; tols=(0.99,0.999,0.9999))
     … +29 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Update(~/.claude/jobs/513a65b3/tmp/profile_component_freeze.jl)
  ⎿  Added 1 line
      1  using PauliStrings
      2 +using LinearAlgebra                                             
      3  using Printf
      4  
      5  function MFIM(N,h)

  Ran 1 shell command (ctrl+o to expand)

● Profiled the component-freeze method across tolerances.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/profile_component_freeze.jl

  Results

  N=16, M=10, steps=20

  ┌────────────────────┬────────┬─────────┬──────────┬───────────┬─────────┐
  │       Method       │  Time  │  Alloc  │  Final   │  rel vs   │  max    │
  │                    │        │         │   len    │    RK4    │  diff   │
  ├────────────────────┼────────┼─────────┼──────────┼───────────┼─────────┤
  │ RK4                │ 0.133s │  0.065  │     1024 │         — │       — │
  │                    │        │      GB │          │           │         │
  ├────────────────────┼────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS tol=0.99 │ 0.552s │  1.400  │     1024 │   1.31e-2 │ 2.24e-1 │
  │                    │        │      GB │          │           │         │
  ├────────────────────┼────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS          │ 1.299s │  2.091  │     1024 │   4.91e-3 │ 8.24e-2 │
  │ tol=0.999          │        │      GB │          │           │         │
  ├────────────────────┼────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS          │ 1.568s │  2.425  │     1024 │   3.41e-3 │ 5.28e-2 │
  │ tol=0.9999         │        │      GB │          │           │         │
  └────────────────────┴────────┴─────────┴──────────┴───────────┴─────────┘

  For smaller N, component-freeze is still slower than RK4.

  N=32, M=12, steps=10

  ┌───────────────────┬─────────┬─────────┬──────────┬───────────┬─────────┐
  │      Method       │  Time   │  Alloc  │  Final   │  rel vs   │  max    │
  │                   │         │         │   len    │    RK4    │  diff   │
  ├───────────────────┼─────────┼─────────┼──────────┼───────────┼─────────┤
  │ RK4               │  0.260s │  0.106  │     4096 │         — │       — │
  │                   │         │      GB │          │           │         │
  ├───────────────────┼─────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS         │ 0.0097s │  0.013  │       95 │   3.44e-3 │ 1.11e-1 │
  │ tol=0.99          │         │      GB │          │           │         │
  ├───────────────────┼─────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS         │ 0.0237s │  0.032  │      235 │   3.96e-3 │ 1.15e-1 │
  │ tol=0.999         │         │      GB │          │           │         │
  ├───────────────────┼─────────┼─────────┼──────────┼───────────┼─────────┤
  │ TrotterTS         │ 0.0585s │  0.094  │      472 │   2.85e-3 │ 7.90e-2 │
  │ tol=0.9999        │         │      GB │          │           │         │
  └───────────────────┴─────────┴─────────┴──────────┴───────────┴─────────┘

  For N=32, component-freeze is much faster and lower allocation, but final
  operator length is drastically smaller than RK4/truncation cap.

  Interpretation

  Component-freeze is acting as a very strong integrated truncation:

  RK4 final len:       4096
  TrotterTS final len: 95–472

  Yet the simple observable trajectory is still within ~0.3–0.4% relative to RK4
  in these short tests.

  This is promising, but we should not claim it is purely a faster exact
  Trotter. It is:

  TS orbit Trotter + component-weight truncation/freezing

  The speedup comes largely from aggressively keeping only high-weight orbit
  components.

  What to check next

  We need to validate this more carefully:

  1. Compare against dense exact for small N, not just RK4.
  2. Sweep componenttol and M.
  3. Track norm/trace_product loss.
  4. Compare more observables than just Xtot.
  5. Check whether freezing is better than dropping.

  The current result suggests component-level truncation is a viable route, but
  the final state is much smaller, so accuracy validation becomes the main
  question.

✻ Sautéed for 2m 9s

❯ Yes perform all these validations                                             

● Write(~/.claude/jobs/513a65b3/tmp/validate_component_freeze.jl)
  ⎿  Wrote 164 lines to
     ../../../.claude/jobs/513a65b3/tmp/validate_component_freeze.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +183 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/validate_component_freeze_corrected.jl)
  ⎿  Wrote 84 lines to ../../../.claude/jobs/513a65b3/tmp/validate_componen
     t_freeze_corrected.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +97 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Performed the requested validation.

  Corrected exact comparison

  The first exact comparison script had a normalization mistake for dense
  observables. I fixed it and reran.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/validate_component_freeze_corrected.jl

  Observables checked:

  Xtot, Ztot, ZZtot, YZtot

  Small systems vs dense exact

  N=6, M=10, T=1.0, dt=0.1

  ┌───────────────────┬───────────┬───────────┬────────────┬───────────┐
  │      Method       │ rel exact │ max exact │ norm ratio │ final len │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ RK4               │   4.44e-4 │   3.14e-3 │   0.999426 │       699 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.99   │   1.63e-2 │   7.60e-2 │        1.0 │        95 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.999  │   5.49e-3 │   2.55e-2 │        1.0 │       233 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.9999 │   5.34e-3 │   2.61e-2 │        1.0 │       419 │
  └───────────────────┴───────────┴───────────┴────────────┴───────────┘

  N=8, M=10, T=1.0, dt=0.1

  ┌───────────────────┬───────────┬───────────┬────────────┬───────────┐
  │      Method       │ rel exact │ max exact │ norm ratio │ final len │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ RK4               │   4.44e-4 │   4.18e-3 │   0.999426 │      1024 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.99   │   1.59e-2 │   9.41e-2 │        1.0 │       114 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.999  │   6.59e-3 │   3.60e-2 │        1.0 │       239 │
  ├───────────────────┼───────────┼───────────┼────────────┼───────────┤
  │ freeze tol=0.9999 │   5.32e-3 │   3.44e-2 │        1.0 │       491 │
  └───────────────────┴───────────┴───────────┴────────────┴───────────┘

  Medium system vs RK4

  N=32, M=12, T=1.0, dt=0.1

  RK4:

  norm_ratio = 0.999426
  final_len = 4096

  ┌──────────────┬──────────┬─────────┬─────────┬───────┬────────┬────────┐
  │    Method    │   time   │  alloc  │ rel vs  │ max   │ norm   │ final  │
  │              │          │         │   RK4   │ diff  │ ratio  │  len   │
  ├──────────────┼──────────┼─────────┼─────────┼───────┼────────┼────────┤
  │ freeze       │ 0.00575s │ 0.207GB │ 2.44e-2 │ 0.652 │    1.0 │     95 │
  │ tol=0.99     │          │         │         │       │        │        │
  ├──────────────┼──────────┼─────────┼─────────┼───────┼────────┼────────┤
  │ freeze       │  0.0124s │ 0.041GB │ 6.68e-3 │ 0.142 │    1.0 │    235 │
  │ tol=0.999    │          │         │         │       │        │        │
  ├──────────────┼──────────┼─────────┼─────────┼───────┼────────┼────────┤
  │ freeze       │  0.0303s │ 0.110GB │ 5.37e-3 │ 0.138 │    1.0 │    472 │
  │ tol=0.9999   │          │         │         │       │        │        │
  └──────────────┴──────────┴─────────┴─────────┴───────┴────────┴────────┘

  Interpretation

  Component-freeze is very fast and aggressively compact, but accuracy is worse
  than RK4 in small exact tests.

  Key observations:

  1. Freezing preserves norm very well:

  norm_ratio = 1.0

  1. That makes sense: low-weight components are carried forward, not dropped.
  2. It produces much smaller operators:

  N=32 final len = 95–472

  2. vs RK4’s 4096.
  3. Accuracy against exact is about:

  ~5e-3 for tol=0.999/0.9999

  3. while RK4 is:

  ~4e-4
  4. Increasing tolerance from 0.999 to 0.9999 does not significantly improve
  small-system exact error. This suggests error is not solely from freezing
  low-weight components; it may also come from:
    - Strang splitting over orbits at dt=0.1,
    - component truncation/freeze interactions,
    - post-step trim,
    - or limitations of this observable set.

  Conclusion

  Component-freeze is a viable aggressive TS-native approximate evolution
  strategy, but not yet an accuracy-equivalent replacement for RK4 at the same
  dt.

  It may still be useful if we compare at matched accuracy/runtime, or if larger
  systems benefit from smaller representations. But for the original goal
  “outperform RK4 in TS case,” we need either:

  - smaller dt / error-controlled component tolerance,
  - better component update for frozen components,
  - or a more accurate low-cost treatment of discarded components.

✻ Cooked for 3m 45s

※ recap: We’re building a genuinely TS-native second-order Trotter method using orbit-level Pauli propagation. Next, test whether component-freezing can be made accurate enough, or replace freezing with a better low-cost update for skipped components.

❯ actually, check if dropping the components and renormalizing is better.       
  frozen components might add more error than missing components                

● That’s a good hypothesis.

  Freezing low-weight components preserves norm, but it also evolves part of the
  operator under the wrong dynamics:

  important components: evolved by orbit flow
  low-weight components: held fixed

  This can introduce phase/dynamical errors that accumulate coherently. Dropping
  them instead gives:

  low-weight components: removed
  then optionally renormalize

  That may be better if the low-weight components are mostly truncation noise.

  We already did a small scratch check for N=6,8 with drop without 
  renormalization, and it was worse than freeze:

  Previous scratch result, N=6

  freeze tol=0.990: rel_rk = 1.64e-2, norm_ratio = 1.00
  drop   tol=0.990: rel_rk = 2.59e-1, norm_ratio = 0.857

  freeze tol=0.999: rel_rk = 5.48e-3, norm_ratio = 1.00
  drop   tol=0.999: rel_rk = 5.91e-2, norm_ratio = 0.982

  Previous scratch result, N=8

  freeze tol=0.990: rel_rk = 1.61e-2, norm_ratio = 1.00
  drop   tol=0.990: rel_rk = 4.86e-1, norm_ratio = 0.885

  freeze tol=0.999: rel_rk = 6.63e-3, norm_ratio = 1.00
  drop   tol=0.999: rel_rk = 8.08e-2, norm_ratio = 0.980

  But that was drop without renormalization. Your suggestion is specifically:

  drop components, then renormalize

  That is worth checking.

  What to test

  For each orbit flow:

  1. Sort components by weight.
  2. Evolve top components.
  3. Drop low-weight components.
  4. Renormalize the result to preserve the pre-flow norm.

  Possible renormalization choices:

  A. Global Hilbert-Schmidt renormalization per orbit flow

  O_new *= sqrt(trace_product(O_old', O_old) / trace_product(O_new', O_new))

  This is simple and likely best for first test.

  B. Coefficient-norm renormalization

  Cheaper but less faithful to TS orbit multiplicities.

  C. No renormalization

  Already tested; poor.

  Expected behavior

  Drop+renorm may improve over drop because it fixes norm loss, but it may still
  be worse than freeze if discarded components carry observable-relevant
  directions.

  The real question:

  Is wrong frozen dynamics worse than missing-but-renormalized components?

  We need data.

  I would test

  For N=6,8 against dense exact and N=32 against RK4:

  freeze
  drop
  drop+renorm

  at:

  tol = 0.99, 0.999, 0.9999

  and measure:

  rel error
  max observable error
  norm ratio
  final length
  runtime

  If drop+renorm wins, we should make the component truncation mode
  configurable:

  componentmode = :freeze | :drop | :renorm

  but only after validating.
✻ Churned for 37s

❯ go ahead. do not wait for verification                                        

● Write(~/.claude/jobs/513a65b3/tmp/compare_component_modes.jl)
  ⎿  Wrote 129 lines to
     ../../../.claude/jobs/513a65b3/tmp/compare_component_modes.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +148 lines (ctrl+o to expand)

● I’ll run a scratch comparison of three component truncation modes now: freeze,
  drop, and drop+renorm, against dense exact for small systems and RK4 for a
  medium case.

  Ran 1 shell command (ctrl+o to expand)

● Tested freeze, drop, and drop+renorm.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/compare_component_modes.jl

  Small systems vs dense exact

  N=6, M=10

  ┌────────┬────────┬───────────┬───────────┬────────────┬─────┬──────────┐
  │  Mode  │  tol   │ rel exact │ max exact │ norm ratio │ len │   time   │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ freeze │   0.99 │   1.63e-2 │   7.60e-2 │        1.0 │  95 │    1.27s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ freeze │  0.999 │   5.49e-3 │   2.55e-2 │        1.0 │ 233 │  0.0095s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ freeze │ 0.9999 │   5.34e-3 │   2.61e-2 │        1.0 │ 419 │   0.019s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ drop   │   0.99 │   2.59e-1 │      1.60 │      0.857 │  12 │ 0.00074s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ drop   │  0.999 │   5.91e-2 │     0.299 │      0.982 │  47 │  0.0022s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ drop   │ 0.9999 │   1.77e-2 │   7.94e-2 │      0.997 │ 116 │  0.0048s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ renorm │   0.99 │   2.54e-1 │      1.60 │        1.0 │  12 │ 0.00079s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ renorm │  0.999 │   5.87e-2 │     0.298 │        1.0 │  47 │  0.0024s │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┼──────────┤
  │ renorm │ 0.9999 │   1.75e-2 │   7.94e-2 │        1.0 │ 116 │  0.0054s │
  └────────┴────────┴───────────┴───────────┴────────────┴─────┴──────────┘

  N=8, M=10

  ┌────────┬────────┬───────────┬───────────┬────────────┬─────┐
  │  Mode  │  tol   │ rel exact │ max exact │ norm ratio │ len │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ freeze │   0.99 │   1.59e-2 │   9.41e-2 │        1.0 │ 114 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ freeze │  0.999 │   6.59e-3 │   3.60e-2 │        1.0 │ 239 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ freeze │ 0.9999 │   5.32e-3 │   3.44e-2 │        1.0 │ 491 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ drop   │   0.99 │   4.86e-1 │      3.33 │      0.885 │   8 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ drop   │  0.999 │   8.08e-2 │     0.412 │      0.980 │  57 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ drop   │ 0.9999 │   5.47e-2 │     0.412 │      0.997 │ 156 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ renorm │   0.99 │   4.85e-1 │      3.33 │        1.0 │   8 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ renorm │  0.999 │   7.88e-2 │     0.412 │        1.0 │  57 │
  ├────────┼────────┼───────────┼───────────┼────────────┼─────┤
  │ renorm │ 0.9999 │   5.47e-2 │     0.412 │        1.0 │ 156 │
  └────────┴────────┴───────────┴───────────┴────────────┴─────┘

  Medium case vs RK4

  N=32, M=12

  ┌────────┬───────┬────────────┬──────────┬────────────┬─────┬─────────┐
  │  Mode  │  tol  │ rel vs RK4 │ max diff │ norm ratio │ len │  time   │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ freeze │  0.99 │    2.44e-2 │    0.652 │        1.0 │  95 │   0.47s │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ freeze │ 0.999 │    6.68e-3 │    0.142 │        1.0 │ 235 │  0.011s │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ drop   │  0.99 │    2.59e-1 │     6.98 │      0.902 │  13 │ 0.0011s │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ drop   │ 0.999 │    6.36e-2 │     1.65 │      0.980 │  17 │ 0.0021s │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ renorm │  0.99 │    2.51e-1 │     6.75 │        1.0 │  13 │ 0.0012s │
  ├────────┼───────┼────────────┼──────────┼────────────┼─────┼─────────┤
  │ renorm │ 0.999 │    6.44e-2 │     1.65 │        1.0 │  17 │ 0.0023s │
  └────────┴───────┴────────────┴──────────┴────────────┴─────┴─────────┘

  Conclusion

  Drop+renorm is not better than freeze.

  Renormalization fixes norm loss, but not the directional/observable error from
  deleting low-weight components. Drop and renorm are much faster because they
  keep very few strings, but they are far less accurate.

  Freeze is clearly better accuracy-wise for the same tolerance.

  So if we keep component-level truncation, the best current mode is freeze.

✻ Cogitated for 1m 22s

❯ ok then instead of complete freeze, use a very cheap update rule for those    
  components. also, if the cost is already way smaller, maybe the tolerance can 
  be set higher to get higher accuracy still with speedups                      

● I’ll replace the low-weight “freeze” tail with a cheap first-order component
  update:

  c_C \leftarrow c_C + dt A_C c_C

  for components below the exact-exponential threshold. Since we already
  build/cache A_C, this should be much cheaper than exp(dt*A_C) and more
  accurate than freezing. Then I’ll sweep stricter tolerances like 0.9999,
  0.99999, maybe 0.999999.

● Update(src/evolution.jl)
  ⎿  Added 4 lines, removed 3 lines
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34  Liouvillians. Each orbit flow is applied by exponentiating conn
          ected components of
      35  the Pauli-orbit graph. Components are processed in descending c
          oefficient weight
      36 -until `componenttol` cumulative weight is reached; lower-weight
         - components are                                                
      37 -frozen for that orbit flow. Components larger than `componentma
         -x` fall back to a                                              
      38 -matrix-free Arnoldi/Krylov exponential action with dimension `k
         -rylovdim` when kept.                                           
      36 +until `componenttol` cumulative weight is reached; lower-weight
         + components receive a                                          
      37 +cheap first-order component update for that orbit flow. Compone
         +nts larger than                                                
      38 +`componentmax` fall back to a matrix-free Arnoldi/Krylov expone
         +ntial action with                                              
      39 +dimension `krylovdim` when kept.                               
      40  """
      41  struct TrotterTS <: AbstractEvolutionMethod
      42      order::Int

● Update(src/evolution.jl)
  ⎿  Added 5 lines, removed 5 lines
      477              end
      478              accumulated += weight
      479          else
      480 -            # Freeze low-weight components: carry their active
          - coefficients forward unchanged.                              
      481 -            for (j, q) in enumerate(plan.component)           
      482 -                !iszero(coeffs[j]) || continue                
      483 -                out = out + typeof(O)([q], [coeffs[j]])       
      484 -            end                                               
      480 +            # Cheap tail update: first-order component flow fo
          +r low-weight components.                                      
      481 +            # This keeps component truncation inside the orbit
          + update without fully freezing                                
      482 +            # the discarded tail.                             
      483 +            coeffs2 = plan.capped ? coeffs : coeffs + dt * (pl
          +an.matrix * coeffs)                                           
      484 +            out = out + typeof(O)(plan.component, coeffs2)    
      485          end
      486      end
      487      return truncation(out)

● Also, just to clarify wording: by “higher tolerance” I’ll test values closer
  to 1.0 (stricter retained weight), e.g. 0.99999, because componenttol means
  cumulative retained weight.

  Ran 1 shell command (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/profile_component_tail_update.jl)
  ⎿  Wrote 62 lines to
     ../../../.claude/jobs/513a65b3/tmp/profile_component_tail_update.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +77 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Implemented the cheap first-order tail update and validated it.

  What changed

  Low-weight components are no longer frozen. They now get a cheap first-order
  component update:

  c_C \leftarrow c_C + dt A_C c_C

  High-weight components still get exact component exponentials.

  So the orbit-flow tail treatment is now:

  important components: exact exp(dt A)
  low-weight components: first-order update

  Tests

  test/evolution.jl passes.

  Validation results

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/profile_component_tail_update.jl

  Small exact: N=6, M=10

  ┌────────────┬────────┬────────┬────────┬─────────┬─────┬───────┬────────┐
  │   Method   │  rel   │  rel   │  max   │  norm   │ len │ time  │ alloc  │
  │            │ exact  │  RK4   │ exact  │  ratio  │     │       │        │
  ├────────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ RK4        │ 4.44e- │      — │ 3.14e- │ 0.99942 │ 699 │     — │      — │
  │            │      4 │        │      3 │       6 │     │       │        │
  ├────────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1      │ 5.36e- │ 5.39e- │ 2.70e- │ 1.00062 │ 699 │ 0.120 │ 0.584G │
  │ tol=0.999  │      3 │      3 │      2 │         │     │     s │      B │
  ├────────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1      │ 5.31e- │ 5.33e- │ 2.62e- │ 1.00007 │ 699 │ 0.109 │ 0.313G │
  │ tol=0.9999 │      3 │      3 │      2 │         │     │     s │      B │
  ├────────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1 tol= │ 5.31e- │ 5.32e- │ 2.61e- │ 1.00001 │ 699 │ 0.121 │ 0.315G │
  │ 0.99999    │      3 │      3 │      2 │         │     │     s │      B │
  └────────────┴────────┴────────┴────────┴─────────┴─────┴───────┴────────┘

  Small exact: N=8, M=10

  ┌───────────┬────────┬────────┬────────┬─────────┬─────┬───────┬────────┐
  │  Method   │  rel   │  rel   │  max   │  norm   │ len │ time  │ alloc  │
  │           │ exact  │  RK4   │ exact  │  ratio  │     │       │        │
  ├───────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ RK4       │ 4.44e- │      — │ 4.18e- │ 0.99942 │ 102 │     — │      — │
  │           │      4 │        │      3 │       6 │   4 │       │        │
  ├───────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1     │ 5.32e- │ 5.34e- │ 3.49e- │ 1.00035 │ 102 │ 0.485 │ 1.337G │
  │ tol=0.999 │      3 │      3 │      2 │         │   4 │     s │      B │
  ├───────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1 tol │ 5.31e- │ 5.32e- │ 3.49e- │ 1.00004 │ 102 │ 0.455 │ 1.162G │
  │ =0.9999   │      3 │      3 │      2 │         │   4 │     s │      B │
  ├───────────┼────────┼────────┼────────┼─────────┼─────┼───────┼────────┤
  │ tail1 tol │ 5.31e- │ 5.32e- │ 3.49e- │     1.0 │ 102 │ 0.498 │ 1.173G │
  │ =0.99999  │      3 │      3 │      2 │         │   4 │     s │      B │
  └───────────┴────────┴────────┴────────┴─────────┴─────┴───────┴────────┘

  Medium: N=32, M=12, compared to RK4

  ┌─────────────────┬─────────┬─────────┬───────────┬──────┬───────┬─────────┐
  │     Method      │ rel RK4 │  max    │   norm    │ len  │ time  │  alloc  │
  │                 │         │  diff   │   ratio   │      │       │         │
  ├─────────────────┼─────────┼─────────┼───────────┼──────┼───────┼─────────┤
  │ tail1 tol=0.999 │ 3.84e-3 │   0.105 │   1.00048 │ 4096 │ 3.97s │  9.82GB │
  ├─────────────────┼─────────┼─────────┼───────────┼──────┼───────┼─────────┤
  │ tail1           │ 3.79e-3 │   0.102 │   1.00005 │ 4096 │ 4.55s │ 11.87GB │
  │ tol=0.9999      │         │         │           │      │       │         │
  ├─────────────────┼─────────┼─────────┼───────────┼──────┼───────┼─────────┤
  │ tail1           │ 3.79e-3 │   0.102 │   1.00001 │ 4096 │ 6.15s │ 11.85GB │
  │ tol=0.99999     │         │         │           │      │       │         │
  └─────────────────┴─────────┴─────────┴───────────┴──────┴───────┴─────────┘

  Interpretation

  The cheap first-order tail update improves the conceptual correctness relative
  to freezing, but performance becomes much worse because it updates all
  low-weight components and keeps the operator at full truncation size.

  Compared to freeze:

  - accuracy is slightly better,
  - norm is close to 1,
  - final length returns to the cap,
  - runtime/allocation blow up.

  So the tail update defeats the component-reduction speedup.

  Conclusion

  First-order tail update is not the right compromise.

  We have three modes:

  freeze: fast, compact, moderate error
  drop/renorm: very fast, too inaccurate
  tail1: more accurate, too slow/heavy

  The practical direction is probably freeze with higher tolerance/validation,
  not tail updates, unless we can update the tail much more cheaply without
  retaining all terms.

✻ Baked for 4m 52s

❯ Ok. Let the action we apply on tail components be something we study later.   
  Let us instead look at what weight truncation we need for RK level accuracy,  
  and whether we still retain a speedup at that point                           

● I’ll revert the production tail action back to freeze for now, then sweep
  componenttol toward 1.0 to find what tolerance is needed to approach RK-level
  accuracy and whether any speedup remains.

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 4 lines
      33  (`order=1`) or second-order (`order=2`, Strang) product over th
          ose TS orbit
      34  Liouvillians. Each orbit flow is applied by exponentiating conn
          ected components of
      35  the Pauli-orbit graph. Components are processed in descending c
          oefficient weight
      36 -until `componenttol` cumulative weight is reached; lower-weight
         - components receive a                                          
      37 -cheap first-order component update for that orbit flow. Compone
         -nts larger than                                                
      38 -`componentmax` fall back to a matrix-free Arnoldi/Krylov expone
         -ntial action with                                              
      39 -dimension `krylovdim` when kept.                               
      36 +until `componenttol` cumulative weight is reached; lower-weight
         + components are                                                
      37 +frozen for that orbit flow. Components larger than `componentma
         +x` fall back to a                                              
      38 +matrix-free Arnoldi/Krylov exponential action with dimension `k
         +rylovdim` when kept.                                           
      39  """
      40  struct TrotterTS <: AbstractEvolutionMethod
      41      order::Int

● Update(src/evolution.jl)
  ⎿  Added 5 lines, removed 5 lines
      476              end
      477              accumulated += weight
      478          else
      479 -            # Cheap tail update: first-order component flow fo
          -r low-weight components.                                      
      480 -            # This keeps component truncation inside the orbit
          - update without fully freezing                                
      481 -            # the discarded tail.                             
      482 -            coeffs2 = plan.capped ? coeffs : coeffs + dt * (pl
          -an.matrix * coeffs)                                           
      483 -            out = out + typeof(O)(plan.component, coeffs2)    
      479 +            # Freeze low-weight components: carry their active
          + coefficients forward unchanged.                              
      480 +            for (j, q) in enumerate(plan.component)           
      481 +                !iszero(coeffs[j]) || continue                
      482 +                out = out + typeof(O)([q], [coeffs[j]])       
      483 +            end                                               
      484          end
      485      end
      486      return truncation(out)

● Write(~/.claude/jobs/513a65b3/tmp/sweep_component_tol_accuracy.jl)
  ⎿  Wrote 40 lines to
     ../../../.claude/jobs/513a65b3/tmp/sweep_component_tol_accuracy.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +49 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Reverted the tail action back to freeze and swept componenttol to see what
  tolerance is needed for RK-level accuracy and whether speedup remains.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/sweep_component_tol_accuracy.jl

  Results

  N=16, M=10, T=1, dt=0.1

  RK4:

  time = 0.0228s
  alloc = 0.027GB
  norm_ratio = 0.999426
  len = 1024

  ┌──────────┬─────────┬────────────┬─────────┬──────────┬────────┬────────┐
  │   tol    │  time   │ speedup vs │  alloc  │ rel vs   │  max   │ final  │
  │          │         │     RK4    │         │   RK4    │  diff  │  len   │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │     0.99 │ 0.0051s │      4.51x │ 0.532GB │  2.44e-2 │  0.326 │     95 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │    0.999 │ 0.0121s │      1.88x │ 0.039GB │  6.68e-3 │ 0.0709 │    235 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │   0.9999 │ 0.0336s │      0.68x │ 0.098GB │  5.37e-3 │ 0.0691 │    472 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │  0.99999 │ 0.0737s │      0.31x │ 0.212GB │  5.34e-3 │ 0.0692 │    797 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │ 0.999999 │  0.230s │      0.10x │ 0.561GB │  5.32e-3 │ 0.0697 │   1024 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │      1.0 │  0.712s │      0.03x │ 1.661GB │  5.32e-3 │ 0.0697 │   1024 │
  └──────────┴─────────┴────────────┴─────────┴──────────┴────────┴────────┘

  N=32, M=12, T=1, dt=0.1

  RK4:

  time = 0.1846s
  alloc = 0.119GB
  norm_ratio = 0.999426
  len = 4096

  ┌──────────┬─────────┬────────────┬─────────┬──────────┬────────┬────────┐
  │   tol    │  time   │ speedup vs │  alloc  │ rel vs   │  max   │ final  │
  │          │         │     RK4    │         │   RK4    │  diff  │  len   │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │     0.99 │ 0.0145s │     12.75x │ 0.207GB │  2.44e-2 │  0.652 │     95 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │    0.999 │ 0.0320s │      5.77x │ 0.044GB │  6.68e-3 │  0.142 │    235 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │   0.9999 │ 0.1057s │      1.75x │ 0.113GB │  5.37e-3 │  0.138 │    472 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │  0.99999 │ 0.1770s │      1.04x │ 0.252GB │  5.34e-3 │  0.138 │    797 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │ 0.999999 │  0.427s │      0.43x │ 0.763GB │  5.32e-3 │  0.139 │   1513 │
  ├──────────┼─────────┼────────────┼─────────┼──────────┼────────┼────────┤
  │      1.0 │   9.84s │      0.02x │  27.2GB │  5.32e-3 │  0.139 │   4096 │
  └──────────┴─────────┴────────────┴─────────┴──────────┴────────┴────────┘

  Interpretation

  There is a speed/accuracy tradeoff, but the error plateaus around:

  rel vs RK4 ≈ 5.3e-3

  even as componenttol -> 1.0.

  So the main remaining error is not from freezing/truncating components.

  It is likely from the orbit-level Strang splitting itself at dt=0.1:

  exp(dt/2 L1) exp(dt/2 L2) ... exp(dt Lm) ...

  compared to RK4’s approximation to the full Hamiltonian flow.

  Increasing componenttol beyond 0.999 mostly increases cost and final length
  without meaningfully improving accuracy.

  Practical operating point

  For N=32:

  tol=0.999:
    speedup ~5.8x over RK4
    rel error vs RK4 ~6.7e-3
    final len 235 vs RK4 4096

  For more accuracy:

  tol=0.9999:
    speedup ~1.75x
    rel error vs RK4 ~5.4e-3
    final len 472

  But going above 0.9999 is not worthwhile at this dt.

  Next question

  If we want RK-level accuracy (~1e-4 vs exact in small tests), raising
  component tolerance will not get us there. We need to reduce Trotter splitting
  error, probably by using smaller dt or a higher-order orbit product formula.

  The component truncation is not the dominant accuracy limiter once tol >= 
  0.999.

✻ Brewed for 3m 18s

❯ I see. What is the dt in rk4? we should use the same                          
● In our comparisons, RK4 and TrotterTS used the same external dt = 0.1.

  But the methods use that dt differently:

  RK4

  For each interval of length dt = 0.1, RK4 performs one fourth-order
  Runge–Kutta step:

  O(t + dt) =
  O + \frac{dt}{6}(k_1 + 2k_2 + 2k_3 + k_4)

  with four full-Hamiltonian commutator evaluations.

  So RK4’s local error is:

  O(dt^5)

  and global error is:

  O(dt^4).

  TrotterTS

  For each interval of length dt = 0.1, TrotterTS performs one second-order
  Strang step over orbit Liouvillians:

  e^{dt/2 L_1}
  e^{dt/2 L_2}
  ...
  e^{dt L_m}
  ...
  e^{dt/2 L_1}.

  So its local splitting error is:

  O(dt^3)

  and global error is:

  O(dt^2).

  Thus using the same dt=0.1 is fair in terms of same step size, but not fair in
  terms of expected accuracy order.

  If we want matched accuracy, TrotterTS likely needs a smaller timestep than
  RK4.

  What we should test

  Use a finer internal Trotter timestep while saving at the same output times.

  Example:

  RK4 save dt = 0.1, one RK4 step per save interval.
  TrotterTS internal dt = 0.05 or 0.025, save every 0.1.

  Since the evolve API uses every tspan interval as an integration step, we can
  just pass a finer times grid to TrotterTS and subsample the history for
  comparison.

  We should compare:

  RK4 dt=0.1
  TrotterTS dt=0.1
  TrotterTS dt=0.05
  TrotterTS dt=0.025

  at matched final/save times, with componenttol fixed, e.g. 0.999 or 0.9999.

  That will tell us whether the observed ~5e-3 error plateau is really Strang
  splitting error.

✻ Brewed for 22s

❯ ok try. my main concern is that we may lose the speedup we gained if we       
  sample more.                                                                  

● Write(~/.claude/jobs/513a65b3/tmp/sweep_trotterts_internal_dt.jl)
  ⎿  Wrote 53 lines to
     ../../../.claude/jobs/513a65b3/tmp/sweep_trotterts_internal_dt.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +58 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Tested whether smaller internal TrotterTS timestep removes the accuracy gap,
  and whether speedup survives.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/sweep_trotterts_internal_dt.jl

  Comparison: RK4 uses save/internal dt=0.1. TrotterTS uses smaller internal dt
  and is subsampled at save times.

  N=16, M=10, componenttol=0.999

  RK4:

  dt=0.1 time=0.0379s len=1024

  ┌─────────────┬────────┬─────────┬────────────┬──────────┬─────────┬─────┐
  │ TS internal │ factor │  time   │ speedup vs │ rel vs   │  max    │ len │
  │      dt     │        │         │     RK4    │   RK4    │  diff   │     │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │         0.1 │      1 │ 0.0147s │      2.57x │  6.68e-3 │ 7.09e-2 │ 235 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │        0.05 │      2 │ 0.0299s │      1.27x │  2.43e-3 │ 2.97e-2 │ 235 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │       0.025 │      4 │ 0.0727s │      0.52x │  1.96e-3 │ 3.25e-2 │ 237 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │      0.0125 │      8 │  0.147s │      0.26x │  2.29e-3 │ 3.25e-2 │ 238 │
  └─────────────┴────────┴─────────┴────────────┴──────────┴─────────┴─────┘

  N=32, M=12, componenttol=0.999

  RK4:

  dt=0.1 time=0.135s len=4096

  ┌─────────────┬────────┬─────────┬────────────┬──────────┬─────────┬─────┐
  │ TS internal │ factor │  time   │ speedup vs │ rel vs   │  max    │ len │
  │      dt     │        │         │     RK4    │   RK4    │  diff   │     │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │         0.1 │      1 │ 0.0136s │      9.90x │  6.68e-3 │ 1.42e-1 │ 235 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │        0.05 │      2 │ 0.0263s │      5.12x │  2.43e-3 │ 5.94e-2 │ 235 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │       0.025 │      4 │ 0.0615s │      2.19x │  1.91e-3 │ 5.96e-2 │ 237 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │      0.0125 │      8 │  0.137s │      0.98x │  2.27e-3 │ 6.09e-2 │ 238 │
  └─────────────┴────────┴─────────┴────────────┴──────────┴─────────┴─────┘

  N=32, M=12, componenttol=0.9999

  RK4:

  dt=0.1 time=0.192s len=4096

  ┌─────────────┬────────┬─────────┬────────────┬──────────┬─────────┬─────┐
  │ TS internal │ factor │  time   │ speedup vs │ rel vs   │  max    │ len │
  │      dt     │        │         │     RK4    │   RK4    │  diff   │     │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │         0.1 │      1 │ 0.0623s │      3.08x │  5.37e-3 │ 1.38e-1 │ 472 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │        0.05 │      2 │ 0.0920s │      2.08x │  1.51e-3 │ 3.47e-2 │ 505 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │       0.025 │      4 │  0.206s │      0.93x │  8.84e-4 │ 2.42e-2 │ 507 │
  ├─────────────┼────────┼─────────┼────────────┼──────────┼─────────┼─────┤
  │      0.0125 │      8 │  0.377s │      0.51x │  8.47e-4 │ 2.47e-2 │ 507 │
  └─────────────┴────────┴─────────┴────────────┴──────────┴─────────┴─────┘

  Interpretation

  Smaller internal timestep does reduce the Strang splitting error, especially
  going from:

  dt=0.1 -> dt=0.05

  At N=32, tol=0.9999, using:

  TrotterTS internal dt = 0.05

  gives:

  ~2.1x speedup over RK4
  rel difference vs RK4 ~1.5e-3
  max diff ~3.5e-2

  At:

  dt = 0.025

  accuracy improves to:

  ~8.8e-4 rel vs RK4

  but speedup is gone.

  So your concern is correct: matching RK-level accuracy eats most or all of the
  speedup.

  Practical candidate

  For a speed/accuracy compromise:

  componenttol = 0.9999
  internal dt = 0.05

  seems promising:

  ~2x faster than RK4
  ~1.5e-3 relative observable difference
  operator length ~505 vs RK4 4096

  For stricter accuracy near 1e-3, speedup is marginal or lost.

✻ Cogitated for 3m 57s

❯  Ok lets revisit. Do we even expect a speedup at rk4 accuracies?              

● Probably not generically at the same accuracy target — at least not from the
  current orbit-Strang + component truncation route.

  The data suggests we can get speedups when we accept a coarser approximation,
  but as we push toward RK4-level accuracy, the speedup mostly disappears.

  Why

  RK4 is fourth order in timestep:

  \text{global error} \sim O(dt^4).

  Our TS orbit Trotter is second order:

  \text{global splitting error} \sim O(dt^2).

  So to reach RK4-level accuracy at dt = 0.1, the Trotter method needs a smaller
  internal timestep.

  That immediately eats speedup.

  What the sweep showed

  For N=32, M=12, componenttol=0.9999:

  ┌───────────────────────┬────────────────┬─────────────────┐
  │ TrotterTS internal dt │ speedup vs RK4 │ rel diff vs RK4 │
  ├───────────────────────┼────────────────┼─────────────────┤
  │                   0.1 │           3.1x │         5.37e-3 │
  ├───────────────────────┼────────────────┼─────────────────┤
  │                  0.05 │           2.1x │         1.51e-3 │
  ├───────────────────────┼────────────────┼─────────────────┤
  │                 0.025 │          0.93x │         8.84e-4 │
  ├───────────────────────┼────────────────┼─────────────────┤
  │                0.0125 │          0.51x │         8.47e-4 │
  └───────────────────────┴────────────────┴─────────────────┘

  So at around 1e-3 relative agreement, speedup is already gone or marginal.

  For tol=0.999, speed is higher, but error plateaus around ~2e-3.

  Why non-TS Trotter can beat RK4

  In the non-TS case, Trotter applies exact single-Pauli rotations. Those
  rotations are very cheap and preserve the exact orbit under a single Pauli
  generator.

  So even with second-order splitting, each substep is extremely cheap.

  In TS orbit Trotter, the orbit generator is:

  H_a^{orb} = \sum_r T_r P_a T_r^\dagger.

  The exact orbit exponential is not cheap unless we exploit more structure. Our
  component method is exact per orbit component, but it creates many components
  and uses truncation/freezing for speed.

  So the cost/accuracy tradeoff is less favorable.

  What we can reasonably expect

  Generic TS only

  If we exploit only translation symmetry, I would not expect a robust speedup
  over RK4 at RK4-level accuracy.

  Translation symmetry gives representation compression, but RK4 already
  benefits from the existing TS commutator. So we are competing against an
  already TS-native RK4.

  The remaining advantage of Trotter would need to come from cheaper finite
  orbit exponentials. But generic orbit exponentials are not cheap enough yet.

  TS + additional structure

  If we allow more structure, then yes, speedups are plausible:

  - commuting orbit terms,
  - colorable local orbits,
  - repeated component signatures,
  - low-rank/closed-form component families,
  - fixed basis/cached component plans across long evolutions.

  But that goes beyond “translation symmetry alone.”

  Honest conclusion

  The original intuition:

  non-TS Trotter beats RK4, so TS Trotter should too

  does not automatically hold, because TS RK4 already uses the compressed TS
  algebra, while current/naive TS Trotter loses the single-Pauli rotation
  advantage.

  A generic TS-native Trotter is mathematically possible, but at RK4 accuracy it
  may not outperform RK4 unless we exploit more than translation symmetry or
  find a much better orbit exponential backend.

  What this means for the project

  The useful outcome may be:

  1. explain why resum-then-Trotter is slow,
  2. formulate the correct TS orbit Trotter,
  3. demonstrate speed/accuracy tradeoffs,
  4. conclude that generic TS-only Trotter is not guaranteed to beat RK4 at RK4
  accuracy,
  5. identify extra conditions under which it can.

  This is a valid research result, even if it means the desired performance goal
  needs refinement.

✻ Cogitated for 34s

※ recap: We’re trying to build a TS-native second-order Trotter that beats RK4 by using orbit-level Pauli propagation. Next, decide whether to continue optimizing component orbit flows or document that generic TS-only Trotter likely won’t beat RK4 at RK4 accuracy.

❯ ok lets clean up the code. Drop the krylov subroutine backend, since we never 
  hit it anyway. make sure that the footprint is as minimal as possible.        
  ⎿  Interrupted · What should Claude do instead?

❯ ok lets clean up the code. Drop the krylov subroutine backend, since we never 
  hit it anyway. make sure that the footprint is as minimal as possible.        

  Read 3 files (ctrl+o to expand)

● I’ll remove the unused Krylov fallback/backend now, since componentmax is high
  enough in our intended use and we don’t want the extra footprint. I’ll also
  simplify the TrotterTS type/docs accordingly and rerun tests.

  Read 1 file (ctrl+o to expand)

● Update(src/evolution.jl)
  ⎿  Added 4 lines, removed 6 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, krylovdim=4, componentmax=1024, compon
         -enttol=0.9999)                                                 
      29 +    TrotterTS(; order=2, componentmax=1024, componenttol=0.9999
         +)                                                              
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
     ...
      34  Liouvillians. Each orbit flow is applied by exponentiating conn
          ected components of
      35  the Pauli-orbit graph. Components are processed in descending c
          oefficient weight
      36  until `componenttol` cumulative weight is reached; lower-weight
           components are
      37 -frozen for that orbit flow. Components larger than `componentma
         -x` fall back to a                                              
      38 -matrix-free Arnoldi/Krylov exponential action with dimension `k
         -rylovdim` when kept.                                           
      37 +frozen for that orbit flow.                                    
      38  """
      39  struct TrotterTS <: AbstractEvolutionMethod
      40      order::Int
      42 -    krylovdim::Int                                             
      41      componentmax::Int
      42      componenttol::Float64
      43  end
      46 -TrotterTS(; order::Integer=2, krylovdim::Integer=4, componentma
         -x::Integer=1024, componenttol::Real=0.9999) =                  
      47 -    TrotterTS(Int(order), Int(krylovdim), Int(componentmax), Fl
         -oat64(componenttol))                                           
      44 +TrotterTS(; order::Integer=2, componentmax::Integer=1024, compo
         +nenttol::Real=0.9999) =                                        
      45 +    TrotterTS(Int(order), Int(componentmax), Float64(componentt
         +ol))                                                           
      46  
      47  """
      48      RK4()

● Update(src/evolution.jl)
  ⎿  Removed 50 lines
      298      return EvolutionResult(collect(tspan), history, OperatorTS
           {Ls,Ps}(Or))
      299  end
      300  
      301 -function _orbit_liouvillian(Ha::Operator{<:PauliStringTS}, O::
          -Operator{<:PauliStringTS}, hbar::Real)                        
      302 -    return orbit_liouvillian(Ha, O; hbar=hbar)                
      303 -end                                                           
      304 -                                                              
      305 -function _axpy_combination(basis, coeffs, template)           
      306 -    O = zero(template)                                        
      307 -    for j in eachindex(coeffs)                                
      308 -        iszero(coeffs[j]) || (O = O + basis[j] * coeffs[j])   
      309 -    end                                                       
      310 -    return O                                                  
      311 -end                                                           
      312 -                                                              
      313 -function _orbit_flow_krylov(Ha::Operator{<:PauliStringTS}, O::
          -Operator{<:PauliStringTS}, dt::Real, hbar::Real, krylovdim::In
          -t, truncation)                                                
      314 -    krylovdim > 0 || throw(ArgumentError("krylovdim must be po
          -sitive"))                                                     
      315 -    β = sqrt(real(trace_product(O', O)))                      
      316 -    β == 0 && return O                                        
      317 -                                                              
      318 -    V = typeof(O)[]                                           
      319 -    Vdag = typeof(O)[]                                        
      320 -    Vnorms = ComplexF64[]                                     
      321 -    v1 = O / β                                                
      322 -    push!(V, v1)                                              
      323 -    push!(Vdag, v1')                                          
      324 -    push!(Vnorms, trace_product(Vdag[1], V[1]))               
      325 -    Hk = zeros(ComplexF64, krylovdim + 1, krylovdim)          
      326 -    actual_dim = krylovdim                                    
      327 -                                                              
      328 -    for j in 1:krylovdim                                      
      329 -        w = truncation(_orbit_liouvillian(Ha, V[j], hbar))    
      330 -        for i in 1:j                                          
      331 -            Hk[i, j] = trace_product(Vdag[i], w) / Vnorms[i]  
      332 -            w = truncation(w - V[i] * Hk[i, j])               
      333 -        end                                                   
      334 -        nw = sqrt(real(trace_product(w', w)))                 
      335 -        if nw < 1e-12                                         
      336 -            actual_dim = j                                    
      337 -            break                                             
      338 -        end                                                   
      339 -        Hk[j + 1, j] = nw                                     
      340 -        v = w / nw                                            
      341 -        push!(V, v)                                           
      342 -        push!(Vdag, v')                                       
      343 -        push!(Vnorms, trace_product(Vdag[end], V[end]))       
      344 -    end                                                       
      345 -                                                              
      346 -    F = exp(dt * Hk[1:actual_dim, 1:actual_dim])              
      347 -    coeffs = β * F[:, 1]                                      
      348 -    return truncation(_axpy_combination(V, coeffs, O))        
      349 -end                                                           
      350 -                                                              
      301  struct _OrbitComponentPlan{P}
      302      component::Vector{P}
      303      index::Dict{P,Int}

● Update(src/evolution.jl)
  ⎿  Error editing file

  Read 1 file (ctrl+o to expand)

● I’m removing the Krylov fallback from TrotterTS entirely, so oversized
  components will now raise a clear error telling the caller to increase
  componentmax.

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      368      return E
      369  end
      370  
      371 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          -, hbar::Real, componentmax::Int, krylovdim::Int, componenttol:
          -:Float64, truncation)                                         
      371 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          +, hbar::Real, componentmax::Int, componenttol::Float64, trunca
          +tion)                                                         
      372      componentmax > 0 || throw(ArgumentError("componentmax must
            be positive"))
      373      0 < componenttol <= 1 || throw(ArgumentError("componenttol
            must be in (0, 1]"))
      374      maxlength = 1000

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      378      end
      379  
      380      assigned = Set{eltype(O.strings)}()
      381 -    component_data = Tuple{Float64,Any,Vector{ComplexF64},Bool
          -}[]                                                           
      381 +    component_data = Tuple{Float64,Any,Vector{ComplexF64}}[]  
      382      total_weight = 0.0
      383  
      384      for q0 in O.strings

● Update(src/evolution.jl)
  ⎿  Added 2 lines, removed 1 line
      397              end
      398          end
      399          has_active || continue
      400 +        plan.capped && throw(ArgumentError("orbit component ex
          +ceeded componentmax=$(componentmax); increase componentmax")) 
      401          total_weight += weight
      401 -        push!(component_data, (weight, plan, coeffs, plan.capp
          -ed))                                                          
      402 +        push!(component_data, (weight, plan, coeffs))         
      403      end
      404  
      405      sort!(component_data; by=x -> x[1], rev=true)

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 15 lines
      407      accumulated = 0.0
      408      out = zero(O)
      409  
      410 -    for (weight, plan, coeffs, capped) in component_data      
      410 +    for (weight, plan, coeffs) in component_data              
      411          if accumulated < keep_weight
      412 -            if capped                                         
      413 -                # Fallback only for kept components that excee
          -d the exact-component cap.                                    
      414 -                seed_strings = eltype(O.strings)[]            
      415 -                seed_coeffs = eltype(O.coeffs)[]              
      416 -                for (j, q) in enumerate(plan.component)       
      417 -                    !iszero(coeffs[j]) || continue            
      418 -                    push!(seed_strings, q)                    
      419 -                    push!(seed_coeffs, coeffs[j])             
      420 -                end                                           
      421 -                out = out + _orbit_flow_krylov(Ha, typeof(O)(s
          -eed_strings, seed_coeffs), dt, hbar, krylovdim, truncation)   
      422 -            else                                              
      423 -                coeffs2 = _component_exp(plan, dt) * coeffs   
      424 -                out = out + typeof(O)(plan.component, coeffs2)
      425 -            end                                               
      412 +            coeffs2 = _component_exp(plan, dt) * coeffs       
      413 +            out = out + typeof(O)(plan.component, coeffs2)    
      414              accumulated += weight
      415          else
      416              # Freeze low-weight components: carry their active
            coefficients forward unchanged.

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      424  end
      425  
      426  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
           TrotterTS, hbar::Real, truncation)
      427 -    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          -d.componentmax, method.krylovdim, method.componenttol, truncat
          -ion)                                                          
      427 +    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          +d.componentmax, method.componenttol, truncation)              
      428  end
      429  
      430  function _orbit_terms(H::Operator{<:PauliStringTS})

  Searched for 1 pattern, ran 1 shell command (ctrl+o to expand)

● Cleaned up the Krylov backend/fallback from the production code.

  Removed

  - _orbit_liouvillian
  - _axpy_combination
  - _orbit_flow_krylov
  - all krylovdim fields/constructor args
  - capped-component Krylov fallback

  rg confirms no remaining krylov, Vdag, Hk, or _axpy_combination in
  src/evolution.jl.

  Current TrotterTS

  Now minimal:

  TrotterTS(; order=2, componentmax=1024, componenttol=0.9999)

  Fields:

  order::Int
  componentmax::Int
  componenttol::Float64

  Backend:

  Strang over TS orbit Liouvillians
  → exact component exponentials for retained components
  → freeze tail components

  If a component exceeds componentmax, it now errors clearly:

  orbit component exceeded componentmax=...; increase componentmax

  No hidden fallback.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

✻ Brewed for 5m 8s
❯ why do we error on components>componentmax? Just continue, even if there is a 
  slow down                                                                     

● You’re right — the error was only a safety guard to prevent accidentally
  building a huge component. If the desired behavior is “keep going even if
  slow,” then componentmax should not be a hard cap. I’ll remove the cap/error
  path so components are always closed exactly.

● Update(src/evolution.jl)
  ⎿  Added 3 lines, removed 4 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 -    TrotterTS(; order=2, componentmax=1024, componenttol=0.9999
         -)                                                              
      29 +    TrotterTS(; order=2, componenttol=0.9999)                  
      30  
      31  Translation-symmetric orbit-level product formula. For `H::Oper
          atorTS`, splits the
      32  Hamiltonian into its representative translation orbits and appl
          ies a first-order
     ...
      38  """
      39  struct TrotterTS <: AbstractEvolutionMethod
      40      order::Int
      41 -    componentmax::Int                                          
      41      componenttol::Float64
      42  end
      44 -TrotterTS(; order::Integer=2, componentmax::Integer=1024, compo
         -nenttol::Real=0.9999) =                                        
      45 -    TrotterTS(Int(order), Int(componentmax), Float64(componentt
         -ol))                                                           
      43 +TrotterTS(; order::Integer=2, componenttol::Real=0.9999) =     
      44 +    TrotterTS(Int(order), Float64(componenttol))               
      45  
      46  """
      47      RK4()

● Update(src/evolution.jl)
  ⎿  Removed 1 line
      302      index::Dict{P,Int}
      303      matrix::Matrix{ComplexF64}
      304      exp_cache::Dict{Float64,Matrix{ComplexF64}}
      305 -    capped::Bool                                              
      305  end
      306  
      307  mutable struct _OrbitFlowCache{P,O}

● Update(src/evolution.jl)
  ⎿  Added 2 lines, removed 3 lines
      311  
      312  _OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlow
           Cache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())
      313  
      314 -function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
          -:PauliStringTS, componentmax::Int, maxlength::Int)            
      314 +function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
          +:PauliStringTS, maxlength::Int)                               
      315      component = typeof(seed)[seed]
      316      index = Dict{typeof(seed),Int}(seed => 1)
      317      queue_index = 1
     ...
      321          edges = orbit_edges(Ha, q; maxlength=maxlength)
      322          for r in edges.strings
      323              haskey(index, r) && continue
      324 -            length(component) >= componentmax && return compon
          -ent, index, true                                              
      324              push!(component, r)
      325              index[r] = length(component)
      326          end
      327      end
      329 -    return component, index, false                            
      328 +    return component, index                                   
      329  end
      330  
      331  function _orbit_component_matrix(Ha::Operator{<:PauliStringTS}
           , component, index, hbar::Real, maxlength::Int)

● Update(src/evolution.jl)
  ⎿  Added 4 lines, removed 4 lines
      342      return A
      343  end
      344  
      345 -function _component_plan(cache::_OrbitFlowCache, Ha::Operator{
          -<:PauliStringTS}, seed::PauliStringTS, hbar::Real, componentma
          -x::Int, maxlength::Int)                                       
      345 +function _component_plan(cache::_OrbitFlowCache, Ha::Operator{
          +<:PauliStringTS}, seed::PauliStringTS, hbar::Real, maxlength::
          +Int)                                                          
      346      plan = get(cache.plans, seed, nothing)
      347      plan !== nothing && return plan
      348  
      349 -    component, index, capped = _orbit_component(Ha, seed, comp
          -onentmax, maxlength)                                          
      350 -    A = capped ? zeros(ComplexF64, 0, 0) : _orbit_component_ma
          -trix(Ha, component, index, hbar, maxlength)                   
      351 -    plan = _OrbitComponentPlan(component, index, A, Dict{Float
          -64,Matrix{ComplexF64}}(), capped)                             
      349 +    component, index = _orbit_component(Ha, seed, maxlength)  
      350 +    A = _orbit_component_matrix(Ha, component, index, hbar, ma
          +xlength)                                                      
      351 +    plan = _OrbitComponentPlan(component, index, A, Dict{Float
          +64,Matrix{ComplexF64}}())                                     
      352      for q in component
      353          cache.plans[q] = plan
      354      end

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 2 lines
      365      return E
      366  end
      367  
      368 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          -, hbar::Real, componentmax::Int, componenttol::Float64, trunca
          -tion)                                                         
      369 -    componentmax > 0 || throw(ArgumentError("componentmax must
          - be positive"))                                               
      368 +function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          +O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          +, hbar::Real, componenttol::Float64, truncation)              
      369      0 < componenttol <= 1 || throw(ArgumentError("componenttol
            must be in (0, 1]"))
      370      maxlength = 1000
      371      coeff_lookup = Dict{eltype(O.strings),eltype(O.coeffs)}()

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      379  
      380      for q0 in O.strings
      381          q0 in assigned && continue
      382 -        plan = _component_plan(cache, Ha, q0, hbar, componentm
          -ax, maxlength)                                                
      382 +        plan = _component_plan(cache, Ha, q0, hbar, maxlength)
      383          coeffs = zeros(ComplexF64, length(plan.component))
      384          weight = 0.0
      385          has_active = false

● Update(src/evolution.jl)
  ⎿  Removed 1 line
      393              end
      394          end
      395          has_active || continue
      396 -        plan.capped && throw(ArgumentError("orbit component ex
          -ceeded componentmax=$(componentmax); increase componentmax")) 
      396          total_weight += weight
      397          push!(component_data, (weight, plan, coeffs))
      398      end

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      419  end
      420  
      421  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
           TrotterTS, hbar::Real, truncation)
      422 -    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          -d.componentmax, method.componenttol, truncation)              
      422 +    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          +d.componenttol, truncation)                                   
      423  end
      424  
      425  function _orbit_terms(H::Operator{<:PauliStringTS})

  Searched for 1 pattern, ran 1 shell command (ctrl+o to expand)

● Removed componentmax and all capped-component logic.

  What changed

  TrotterTS is now:

  TrotterTS(; order=2, componenttol=0.9999)

  Fields:

  order::Int
  componenttol::Float64

  No componentmax, no capped, no fallback, no errors for large components.
  Components are closed fully and exponentiated, even if that is slow.

  Verified cleanup

  Searched for:

  componentmax
  capped
  krylov

  No matches remain in src/evolution.jl or test/evolution.jl.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

✻ Cooked for 3m 32s

※ recap: We’re building a generic TS-native Trotter method using orbit-component exponentials with component-weight truncation. Next, verify the cleaned implementation’s accuracy/speed tradeoff after removing Krylov and component caps.

❯ why do we have two _orbit_flow functions? just keep one. same with            
  _active_shift_key                                                             

  Searched for 1 pattern, read 2 files (ctrl+o to expand)

● Agreed. I’ll inline/remove the thin wrappers so there is only one orbit-flow
  helper and no unnecessary active-shift key helper.

● Update(src/evolution.jl)
  ⎿  Added 2 lines, removed 1 line
      365      return E
      366  end
      367  
      368 -function _orbit_flow_component(Ha::Operator{<:PauliStringTS}, 
          -O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache
          -, hbar::Real, componenttol::Float64, truncation)              
      368 +function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
          +r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
          +TrotterTS, hbar::Real, truncation)                            
      369 +    componenttol = method.componenttol                        
      370      0 < componenttol <= 1 || throw(ArgumentError("componenttol
            must be in (0, 1]"))
      371      maxlength = 1000
      372      coeff_lookup = Dict{eltype(O.strings),eltype(O.coeffs)}()

● Update(src/evolution.jl)
  ⎿  Removed 3 lines
      419      return truncation(out)
      420  end
      421  
      422 -function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
          -r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
          -TrotterTS, hbar::Real, truncation)                            
      423 -    return _orbit_flow_component(Ha, O, dt, cache, hbar, metho
          -d.componenttol, truncation)                                   
      424 -end                                                           
      422  
      423  function _orbit_terms(H::Operator{<:PauliStringTS})
      424      terms = typeof(H)[]

● Update(src/translation_symmetry.jl)
  ⎿  Removed 4 lines
      291      return valid, shifts
      292  end
      293  
      294 -function _active_shift_key(shifts::Tuple, Ls::Tuple)          
      295 -    return _linear_site(shifts, Ls)                           
      296 -end                                                           
      297 -                                                              
      294  function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStr
           ingTS; epsilon::Real=0, maxlength::Int=1000)
      295      checklength(A, q)
      296      Ls = qubitsize(A)

● Update(src/translation_symmetry.jl)
  ⎿  Added 1 line, removed 1 line
      318  
      319                  valid, shifts = _active_shift_tuple(site_p, si
           te_q, Ls, Ps)
      320                  valid || continue
      321 -                key = _active_shift_key(shifts, Ls)           
      321 +                key = _linear_site(shifts, Ls)                
      322                  seen[key] && continue
      323                  seen[key] = true
      324  

  Ran 1 shell command (ctrl+o to expand)

● Cleaned both up.

  _orbit_flow

  There is now only one orbit-flow helper:

  function _orbit_flow(...)

  Removed the extra _orbit_flow_component wrapper.

  _active_shift_key

  Removed _active_shift_key entirely and replaced its only use with:

  key = _linear_site(shifts, Ls)

  Verified

  Search now shows only one orbit-flow function and no _active_shift_key:

  src/evolution.jl:368:function _orbit_flow(...)

  Tests pass:

  evolve | 8 passed

✻ Brewed for 2m 13s

❯ consider moving orbit functions into a different file                         
● Write(src/orbit.jl)
  ⎿  Wrote 80 lines to src/orbit.jl
       1 function _coord(site::Integer, Ls::Tuple, dim::Integer)
       2     x = site - 1
       3     for k in 1:(dim - 1)
       4         x = div(x, Ls[k])
       5     end
       6     return mod(x, Ls[dim]) + 1
       7 end
       8 
       9 function _linear_site(coords::Tuple, Ls::Tuple)
      10     site = 1
     … +79 lines (ctrl+o to expand)

● Yes — that will make the scope clearer. I’ll move the generic orbit graph
  helpers (orbit_edges, orbit_liouvillian, and their small coordinate helpers)
  out of translation_symmetry.jl into a new src/orbit.jl, then include it after
  translation_symmetry.jl.

● Update(src/translation_symmetry.jl)
  ⎿  Error editing file

  Read 2 files (ctrl+o to expand)

● Update(src/translation_symmetry.jl)
  ⎿  Removed 90 lines
      257      return op
      258  end
      259  
      260 -@inline function _coord(site::Integer, Ls::Tuple, dim::Integer
          -)                                                             
      261 -    x = site - 1                                              
      262 -    for k in 1:(dim - 1)                                      
      263 -        x = div(x, Ls[k])                                     
      264 -    end                                                       
      265 -    return mod(x, Ls[dim]) + 1                                
      266 -end                                                           
      267 -                                                              
      268 -@inline function _linear_site(coords::Tuple, Ls::Tuple)       
      269 -    site = 1                                                  
      270 -    stride = 1                                                
      271 -    for k in eachindex(Ls)                                    
      272 -        site += (coords[k] - 1) * stride                      
      273 -        stride *= Ls[k]                                       
      274 -    end                                                       
      275 -    return site                                               
      276 -end                                                           
      277 -                                                              
      278 -function _active_shift_tuple(site_p::Integer, site_q::Integer,
          - Ls::Tuple, Ps::Tuple)                                        
      279 -    valid = true                                              
      280 -    shifts = ntuple(length(Ls)) do k                          
      281 -        cp = _coord(site_p, Ls, k)                            
      282 -        cq = _coord(site_q, Ls, k)                            
      283 -        if Ps[k]                                              
      284 -            d = mod(cp - cq, Ls[k])                           
      285 -            iszero(d) ? Ls[k] : d                             
      286 -        else                                                  
      287 -            valid &= cp == cq                                 
      288 -            1                                                 
      289 -        end                                                   
      290 -    end                                                       
      291 -    return valid, shifts                                      
      292 -end                                                           
      293 -                                                              
      294 -function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStr
          -ingTS; epsilon::Real=0, maxlength::Int=1000)                  
      295 -    checklength(A, q)                                         
      296 -    Ls = qubitsize(A)                                         
      297 -    Ps = periodicflags(A)                                     
      298 -    N = Base.prod(Ls)                                         
      299 -    d = emptydict(A)                                          
      300 -    qrep = representative(q)                                  
      301 -    uq = qrep.v | qrep.w                                      
      302 -    iszero(uq) && return typeof(A)(collect(keys(d)), collect(v
          -alues(d)))                                                    
      303 -                                                              
      304 -    for (p, c) in zip(A.strings, A.coeffs)                    
      305 -        prep = representative(p)                              
      306 -        up = prep.v | prep.w                                  
      307 -        iszero(up) && continue                                
      308 -        seen = falses(N)                                      
      309 -        xp = up                                               
      310 -        while !iszero(xp)                                     
      311 -            site_p = trailing_zeros(xp) + 1                   
      312 -            xp &= xp - one(xp)                                
      313 -                                                              
      314 -            xq = uq                                           
      315 -            while !iszero(xq)                                 
      316 -                site_q = trailing_zeros(xq) + 1               
      317 -                xq &= xq - one(xq)                            
      318 -                                                              
      319 -                valid, shifts = _active_shift_tuple(site_p, si
          -te_q, Ls, Ps)                                                 
      320 -                valid || continue                             
      321 -                key = _linear_site(shifts, Ls)                
      322 -                seen[key] && continue                         
      323 -                seen[key] = true                              
      324 -                                                              
      325 -                out, k = commutator(prep, shift(qrep, Ls, Ps, 
          -shifts))                                                      
      326 -                coeff = c * k                                 
      327 -                if (k != 0) && (abs(coeff) > epsilon) && pauli
          -_weight(out) < maxlength                                      
      328 -                    setwith!(+, d, PauliStringTS{Ls,Ps}(out), 
          -coeff)                                                        
      329 -                end                                           
      330 -            end                                               
      331 -        end                                                   
      332 -    end                                                       
      333 -    return typeof(A)(collect(keys(d)), collect(values(d)))    
      334 -end                                                           
      335 -                                                              
      336 -function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Op
          -erator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlen
          -gth::Int=1000)                                                
      337 -    checklength(A, B)                                         
      338 -    d = emptydict(A)                                          
      339 -    for (q, cq) in zip(B.strings, B.coeffs)                   
      340 -        edges = orbit_edges(A, q; epsilon=epsilon, maxlength=m
          -axlength)                                                     
      341 -        for (p, c) in zip(edges.strings, edges.coeffs)        
      342 -            setwith!(+, d, p, cq * c)                         
      343 -        end                                                   
      344 -    end                                                       
      345 -    o = typeof(A)(collect(keys(d)), collect(values(d)))       
      346 -    o = (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) 
          -: o                                                           
      347 -    return 1im * o / hbar                                     
      348 -end                                                           
      349 -                                                              
      260  function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Op
           erator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)
      261      checklength(A, B)
      262      Ls = qubitsize(A)

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line
      39  include("paulistring.jl")
      40  include("operator.jl")
      41  include("translation_symmetry.jl")
      42 +include("orbit.jl")                                            
      43  include("io.jl")
      44  include("operations.jl")
      45  include("operations_strings.jl")

  Ran 1 shell command (ctrl+o to expand)

● Moved the orbit graph helpers into a separate file.

  New file

  src/orbit.jl

  Contains:

  - _coord
  - _linear_site
  - _active_shift_tuple
  - orbit_edges
  - orbit_liouvillian

  Updated include order

  src/PauliStrings.jl now has:

  include("translation_symmetry.jl")
  include("orbit.jl")
  include("io.jl")

  Removed from translation_symmetry.jl

  The orbit helper functions are no longer in translation_symmetry.jl.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

✻ Baked for 2m 56s

❯ what about the other _orbit functions?                                        

  Searched for 1 pattern (ctrl+o to expand)

● Good point. The component/cache helpers are also orbit-graph machinery and
  should live with orbit_edges. I’ll move the generic ones to src/orbit.jl and
  leave only the evolution-method glue (_orbit_flow, _orbit_terms) in
  evolution.jl, since those depend on TrotterTS.

  Read 2 files (ctrl+o to expand)

● Update(src/orbit.jl)
  ⎿  Added 68 lines
       78      checklength(A, B)
       79      return 1im * commutator(A, B; epsilon=epsilon, maxlength=m
           axlength) / hbar
       80  end
       81 +                                                              
       82 +struct _OrbitComponentPlan{P}                                 
       83 +    component::Vector{P}                                      
       84 +    index::Dict{P,Int}                                        
       85 +    matrix::Matrix{ComplexF64}                                
       86 +    exp_cache::Dict{Float64,Matrix{ComplexF64}}               
       87 +end                                                           
       88 +                                                              
       89 +mutable struct _OrbitFlowCache{P,O}                           
       90 +    plans::Dict{P,_OrbitComponentPlan{P}}                     
       91 +    fallback::Dict{P,O}                                       
       92 +end                                                           
       93 +                                                              
       94 +_OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlow
          +Cache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())     
       95 +                                                              
       96 +function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
          +:PauliStringTS, maxlength::Int)                               
       97 +    component = typeof(seed)[seed]                            
       98 +    index = Dict{typeof(seed),Int}(seed => 1)                 
       99 +    queue_index = 1                                           
      100 +    while queue_index <= length(component)                    
      101 +        q = component[queue_index]                            
      102 +        queue_index += 1                                      
      103 +        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      104 +        for r in edges.strings                                
      105 +            haskey(index, r) && continue                      
      106 +            push!(component, r)                               
      107 +            index[r] = length(component)                      
      108 +        end                                                   
      109 +    end                                                       
      110 +    return component, index                                   
      111 +end                                                           
      112 +                                                              
      113 +function _orbit_component_matrix(Ha::Operator{<:PauliStringTS}
          +, component, index, hbar::Real, maxlength::Int)               
      114 +    n = length(component)                                     
      115 +    A = zeros(ComplexF64, n, n)                               
      116 +    for (j, q) in enumerate(component)                        
      117 +        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      118 +        for (r, c) in zip(edges.strings, edges.coeffs)        
      119 +            i = get(index, r, 0)                              
      120 +            i == 0 && throw(ArgumentError("orbit component is 
          +not closed"))                                                 
      121 +            A[i, j] += 1im * c / hbar                         
      122 +        end                                                   
      123 +    end                                                       
      124 +    return A                                                  
      125 +end                                                           
      126 +                                                              
      127 +function _component_plan(cache::_OrbitFlowCache, Ha::Operator{
          +<:PauliStringTS}, seed::PauliStringTS, hbar::Real, maxlength::
          +Int)                                                          
      128 +    plan = get(cache.plans, seed, nothing)                    
      129 +    plan !== nothing && return plan                           
      130 +                                                              
      131 +    component, index = _orbit_component(Ha, seed, maxlength)  
      132 +    A = _orbit_component_matrix(Ha, component, index, hbar, ma
          +xlength)                                                      
      133 +    plan = _OrbitComponentPlan(component, index, A, Dict{Float
          +64,Matrix{ComplexF64}}())                                     
      134 +    for q in component                                        
      135 +        cache.plans[q] = plan                                 
      136 +    end                                                       
      137 +    return plan                                               
      138 +end                                                           
      139 +                                                              
      140 +function _component_exp(plan::_OrbitComponentPlan, dt::Real)  
      141 +    key = Float64(dt)                                         
      142 +    E = get(plan.exp_cache, key, nothing)                     
      143 +    if E === nothing                                          
      144 +        E = exp(key * plan.matrix)                            
      145 +        plan.exp_cache[key] = E                               
      146 +    end                                                       
      147 +    return E                                                  
      148 +end                                                           

● Update(src/evolution.jl)
  ⎿  Removed 68 lines
      297      return EvolutionResult(collect(tspan), history, OperatorTS
           {Ls,Ps}(Or))
      298  end
      299  
      300 -struct _OrbitComponentPlan{P}                                 
      301 -    component::Vector{P}                                      
      302 -    index::Dict{P,Int}                                        
      303 -    matrix::Matrix{ComplexF64}                                
      304 -    exp_cache::Dict{Float64,Matrix{ComplexF64}}               
      305 -end                                                           
      306 -                                                              
      307 -mutable struct _OrbitFlowCache{P,O}                           
      308 -    plans::Dict{P,_OrbitComponentPlan{P}}                     
      309 -    fallback::Dict{P,O}                                       
      310 -end                                                           
      311 -                                                              
      312 -_OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlow
          -Cache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())     
      313 -                                                              
      314 -function _orbit_component(Ha::Operator{<:PauliStringTS}, seed:
          -:PauliStringTS, maxlength::Int)                               
      315 -    component = typeof(seed)[seed]                            
      316 -    index = Dict{typeof(seed),Int}(seed => 1)                 
      317 -    queue_index = 1                                           
      318 -    while queue_index <= length(component)                    
      319 -        q = component[queue_index]                            
      320 -        queue_index += 1                                      
      321 -        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      322 -        for r in edges.strings                                
      323 -            haskey(index, r) && continue                      
      324 -            push!(component, r)                               
      325 -            index[r] = length(component)                      
      326 -        end                                                   
      327 -    end                                                       
      328 -    return component, index                                   
      329 -end                                                           
      330 -                                                              
      331 -function _orbit_component_matrix(Ha::Operator{<:PauliStringTS}
          -, component, index, hbar::Real, maxlength::Int)               
      332 -    n = length(component)                                     
      333 -    A = zeros(ComplexF64, n, n)                               
      334 -    for (j, q) in enumerate(component)                        
      335 -        edges = orbit_edges(Ha, q; maxlength=maxlength)       
      336 -        for (r, c) in zip(edges.strings, edges.coeffs)        
      337 -            i = get(index, r, 0)                              
      338 -            i == 0 && throw(ArgumentError("orbit component is 
          -not closed"))                                                 
      339 -            A[i, j] += 1im * c / hbar                         
      340 -        end                                                   
      341 -    end                                                       
      342 -    return A                                                  
      343 -end                                                           
      344 -                                                              
      345 -function _component_plan(cache::_OrbitFlowCache, Ha::Operator{
          -<:PauliStringTS}, seed::PauliStringTS, hbar::Real, maxlength::
          -Int)                                                          
      346 -    plan = get(cache.plans, seed, nothing)                    
      347 -    plan !== nothing && return plan                           
      348 -                                                              
      349 -    component, index = _orbit_component(Ha, seed, maxlength)  
      350 -    A = _orbit_component_matrix(Ha, component, index, hbar, ma
          -xlength)                                                      
      351 -    plan = _OrbitComponentPlan(component, index, A, Dict{Float
          -64,Matrix{ComplexF64}}())                                     
      352 -    for q in component                                        
      353 -        cache.plans[q] = plan                                 
      354 -    end                                                       
      355 -    return plan                                               
      356 -end                                                           
      357 -                                                              
      358 -function _component_exp(plan::_OrbitComponentPlan, dt::Real)  
      359 -    key = Float64(dt)                                         
      360 -    E = get(plan.exp_cache, key, nothing)                     
      361 -    if E === nothing                                          
      362 -        E = exp(key * plan.matrix)                            
      363 -        plan.exp_cache[key] = E                               
      364 -    end                                                       
      365 -    return E                                                  
      366 -end                                                           
      367 -                                                              
      300  function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operato
           r{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::
           TrotterTS, hbar::Real, truncation)
      301      componenttol = method.componenttol
      302      0 < componenttol <= 1 || throw(ArgumentError("componenttol
            must be in (0, 1]"))

  Ran 1 shell command (ctrl+o to expand)

● Moved the remaining generic orbit component/cache helpers into src/orbit.jl.

  Now in src/orbit.jl

  - orbit_edges
  - orbit_liouvillian
  - _OrbitComponentPlan
  - _OrbitFlowCache
  - _orbit_component
  - _orbit_component_matrix
  - _component_plan
  - _component_exp
  - coordinate/active-shift helpers

  Still in src/evolution.jl

  Only the evolution-specific glue remains:

  - _orbit_flow
  - _orbit_terms
  - _trotterts_step
  - _evolve(::TrotterTS, ...)

  This keeps evolution.jl focused on time stepping and orbit.jl focused on
  Pauli-orbit graph/component machinery.

  Tests

  test/evolution.jl passes:

  evolve | 8 passed

✻ Crunched for 2m 18s

❯ also move _orbit_flow and _orbit_terms?                                       

  Read 1 file (ctrl+o to expand)

● I would move _orbit_terms, yes. But I’d keep _orbit_flow in evolution.jl
  unless we rename/generalize it.

  Reason:

  - _orbit_terms(H) is purely a helper for decomposing an OperatorTS into
  single-orbit terms. It belongs naturally in orbit.jl.
  - _orbit_flow(...) currently takes method::TrotterTS, reads
  method.componenttol, applies the method’s truncation policy, and is tightly
  coupled to the evolution method. It is not just orbit graph machinery; it is
  the TrotterTS policy.

  A cleaner split would be:

  # orbit.jl
  _orbit_terms(H)
  _orbit_component_flow(Ha, O, dt, cache; componenttol, hbar, truncation)

  # evolution.jl
  _orbit_flow(Ha, O, dt, cache, method, hbar, truncation) =
      _orbit_component_flow(Ha, O, dt, cache;
                            componenttol=method.componenttol,
                            hbar=hbar,
                            truncation=truncation)

  That would keep policy in evolution.jl and mechanics in orbit.jl.

  If you want minimal footprint, I can move both, but then orbit.jl will depend
  on TrotterTS, which is defined in evolution.jl and included later — so that
  direction is awkward unless we either:

  1. move TrotterTS type earlier, or
  2. make _orbit_flow not depend on TrotterTS.

  I recommend option 2: move the mechanical function under a name like
  _orbit_component_flow, and leave the tiny method wrapper in evolution.jl.

✻ Baked for 48s

❯ move _orbit_terms but keep _orbit_flow                                        

  Read 1 file (ctrl+o to expand)

● Update(src/orbit.jl)
  ⎿  Added 8 lines
      146      end
      147      return E
      148  end
      149 +                                                              
      150 +function _orbit_terms(H::Operator{<:PauliStringTS})           
      151 +    terms = typeof(H)[]                                       
      152 +    for j in 1:length(H)                                      
      153 +        push!(terms, typeof(H)([H.strings[j]], [H.coeffs[j]]))
      154 +    end                                                       
      155 +    return terms                                              
      156 +end                                                           

● Update(src/evolution.jl)
  ⎿  Removed 9 lines
      352  end
      353  
      354  
      355 -function _orbit_terms(H::Operator{<:PauliStringTS})           
      356 -    terms = typeof(H)[]                                       
      357 -    for j in 1:length(H)                                      
      358 -        c, p = H[j]                                           
      359 -        push!(terms, typeof(H)([p], [H.coeffs[j]]))           
      360 -    end                                                       
      361 -    return terms                                              
      362 -end                                                           
      363 -                                                              
      355  function _trotterts_step(Hterms, caches, O::Operator{<:PauliSt
           ringTS}, dt::Real, method::TrotterTS, hbar::Real, truncation)
      356      method.order ∈ (1, 2) || throw(ArgumentError("order must b
           e 1 or 2, got $(method.order)"))
      357      isempty(Hterms) && return O

  Ran 2 shell commands (ctrl+o to expand)
● Moved _orbit_terms into src/orbit.jl and kept _orbit_flow in src/evolution.jl,
  as requested.

  Current placement:

  src/orbit.jl:
    _orbit_terms
    orbit_edges
    orbit_liouvillian
    _OrbitComponentPlan
    _OrbitFlowCache
    _orbit_component
    _orbit_component_matrix
    _component_plan
    _component_exp

  src/evolution.jl:
    _orbit_flow
    _trotterts_step
    _evolve(::TrotterTS, ...)

  Tests pass:

  evolve | 8 passed

✻ Cooked for 2m 6s

※ recap: We’re building a generic TS-native Trotter method using orbit-component exponentials with component-weight freezing. Next, we should finalize cleanup/tests and decide whether to keep the current component-freeze implementation based on accuracy/speed tradeoffs.

❯ make a commit, with description                                               

  Ran 3 shell commands (ctrl+o to expand)

● Full test suite passes after the cleanup/refactor:

  PauliStrings tests passed

  Current production state:

  - TrotterTS(; order=2, componenttol=0.9999)
  - no Krylov backend/fallback
  - no componentmax
  - orbit graph/component machinery in src/orbit.jl
  - evolution policy in src/evolution.jl
  - component freeze tail behavior retained
  - full tests green

✻ Churned for 2m 11s

※ recap: We’re building a generic TS-native Trotter method using orbit-level component exponentials with component-weight truncation. Next, finalize the cleaned implementation and commit the current tested state.

❯ make a git commit, with description                                           

  Committed beb5526, ran 1 shell command (ctrl+o to expand)

● Committed the refined TS orbit Trotter work.

  beb5526 Refine TS orbit Trotter flow

  Commit includes:

  - TS orbit component evolution with component-weight truncation
  - dedicated src/orbit.jl
  - cleanup of Krylov/component cap paths
  - theory document updates
  - attempts flow document

  Commit message includes the test command:

  julia --project=/home/dhruva/Programs/unitaryHack/PauliStrings.jl -e 'using 
  Pkg; Pkg.test()'

✻ Sautéed for 1m 4s

※ recap: We’re building a TS-native Trotter method using orbit-component evolution with component-weight truncation. Next, validate the committed implementation’s accuracy/speed tradeoff more broadly and decide if the current approximation is acceptable.

❯ Ok, now perform a wide variety of validation checks, and explicitly list RK4  
  time, TrotterTS time, and allocs for both.                                    

● I’ll run a validation sweep across small exact-checkable systems and medium TS
  systems, with RK4 and TrotterTS on the same time grids. I’ll report wall
  time, allocations, final lengths, norm ratios, and observable errors.

● Write(~/.claude/jobs/513a65b3/tmp/validation_sweep_trotterts.jl)
  ⎿  Wrote 94 lines to
     ../../../.claude/jobs/513a65b3/tmp/validation_sweep_trotterts.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function MFIM(N,h)
       6     H=Operator(N)
       7     H += -h,"X",1
       8     H += -h/2,"Z",2
       9     H += "Z",1,"Z",2
      10     OperatorTS{(N,)}(H)
     … +114 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Ran a broader validation sweep comparing RK4 and TrotterTS times/allocations
  explicitly.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/validation_sweep_trotterts.jl

  Observables checked:

  Xtot, Ztot, ZZtot, YZtot, XZtot

  Small exact-checkable cases

  MFIM, N=4, M=8, T=1, dt=0.1, tol=0.9999

  ┌───────────┬──────────┬─────────┬───────────┬────────────┐
  │  Method   │   Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼──────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.00110s │ 0.002GB │        69 │   0.999426 │
  ├───────────┼──────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.00424s │ 0.011GB │        68 │        1.0 │
  └───────────┴──────────┴─────────┴───────────┴────────────┘

  Errors:

  TS vs RK4        rel=5.52e-3, max=1.72e-2
  RK4 vs exact     rel=5.07e-4, max=2.10e-3
  TrotterTS exact  rel=5.46e-3, max=1.72e-2

  MFIM, N=6, M=10, T=1, dt=0.1, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.0100s │ 0.013GB │       699 │   0.999426 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0316s │ 0.080GB │       427 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Errors:

  TS vs RK4        rel=5.41e-3, max=2.59e-2
  RK4 vs exact     rel=5.06e-4, max=3.14e-3
  TrotterTS exact  rel=5.37e-3, max=2.59e-2

  TFIM, N=8, M=10, T=1, dt=0.1, tol=0.9999

  ┌───────────┬───────────┬─────────┬───────────┬────────────┐
  │  Method   │   Time    │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼───────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.000475s │ 0.001GB │        30 │   0.999501 │
  ├───────────┼───────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.000712s │ 0.002GB │        22 │        1.0 │
  └───────────┴───────────┴─────────┴───────────┴────────────┘

  Errors:

  TS vs RK4        rel=5.35e-3, max=3.55e-2
  RK4 vs exact     rel=4.01e-4, max=3.92e-3
  TrotterTS exact  rel=5.31e-3, max=3.57e-2

  MFIM, N=8, M=10, T=1, dt=0.05, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.0438s │ 0.059GB │      1024 │   0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0725s │ 0.193GB │       504 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Errors:

  TS vs RK4        rel=1.61e-3, max=8.42e-3
  RK4 vs exact     rel=3.16e-5, max=2.50e-4
  TrotterTS exact  rel=1.60e-3, max=8.43e-3

  Medium cases vs RK4

  MFIM, N=16, M=10, T=1, dt=0.1, tol=0.999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.0248s │ 0.028GB │      1024 │   0.999426 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0109s │ 0.032GB │       210 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Speedup = 2.27x
  TS vs RK4 rel=5.42e-3, max=7.47e-2

  MFIM, N=16, M=10, T=1, dt=0.05, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │ 0.0539s │ 0.065GB │      1024 │   0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0758s │ 0.222GB │       493 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Speedup = 0.71x
  TS vs RK4 rel=1.39e-3, max=1.71e-2

  MFIM, N=32, M=12, T=1, dt=0.1, tol=0.999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │  0.128s │ 0.123GB │      4096 │   0.999426 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0124s │ 0.038GB │       210 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Speedup = 10.37x
  TS vs RK4 rel=5.42e-3, max=1.50e-1

  MFIM, N=32, M=12, T=1, dt=0.05, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬────────────┐
  │  Method   │  Time   │  Alloc  │ Final len │ Norm ratio │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ RK4       │  0.322s │ 0.298GB │      4096 │   0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼────────────┤
  │ TrotterTS │ 0.0849s │ 0.269GB │       493 │        1.0 │
  └───────────┴─────────┴─────────┴───────────┴────────────┘

  Speedup = 3.80x
  TS vs RK4 rel=1.35e-3, max=3.41e-2

  Summary

  TrotterTS gives real speedups for medium TS cases when accepting ~`1e-3–5e-3`
  observable-level error relative to RK4.

  But for small exact cases, RK4 is much more accurate at the same timestep:

  RK4 exact error ~ 1e-4 to 1e-5
  TrotterTS exact error ~ 1e-3 to 5e-3

  The practical tradeoff looks like:

  - dt=0.1, tol=0.999: large speedups, ~`5e-3` relative observable error.
  - dt=0.05, tol=0.9999: still speedup at N=32, ~`1.3e-3` relative observable
  error.
  - RK4-level exact accuracy is not matched yet.

✻ Brewed for 2m 15s

※ recap: We’re building a TS-native Trotter method using orbit-component exponentials with component-weight truncation. Next, validate the speed/accuracy tradeoff further and decide whether this approximation is acceptable for the PR.

❯ Did we also implement similar component applications?                         

● Yes. The current TrotterTS implementation is based on component applications.

  For each orbit Hamiltonian term Ha, the orbit flow:

  O \mapsto e^{dt L_{Ha}} O

  is applied by decomposing the current OperatorTS into connected components of
  the Pauli-orbit graph, then applying the component action.

  Where it is implemented

  The graph/component machinery is in:

  src/orbit.jl

  Key pieces:

  orbit_edges(Ha, q)

  Builds outgoing graph edges from one TS Pauli orbit q.

  _orbit_component(Ha, seed, maxlength)

  Builds the connected component containing seed.

  _orbit_component_matrix(Ha, component, index, hbar, maxlength)

  Builds the finite Liouvillian matrix for that component.

  _component_plan(cache, Ha, seed, hbar, maxlength)

  Builds/caches the component and its matrix.

  _component_exp(plan, dt)

  Caches and returns:

  exp(dt * plan.matrix)

  How TrotterTS uses it

  The evolution glue is in:

  src/evolution.jl

  In _orbit_flow, for each orbit flow:

  1. Build coefficient lookup for the current O.
  2. Find each active component using _component_plan.
  3. Compute each component’s active coefficient weight.
  4. Sort components by descending weight.
  5. For high-weight components:
    - apply exact component exponential:

  coeffs2 = _component_exp(plan, dt) * coeffs
  6. For low-weight components:
    - freeze them unchanged.

  So yes, the applied update is component-based.

  What we did not implement

  We did not keep:

  - Krylov component fallback,
  - Taylor-2 component update,
  - drop/renormalize tail update,
  - first-order tail update.

  The current production behavior is:

  important components: exact local component exponential
  low-weight components: frozen

  controlled by:

  componenttol

  Default:

  TrotterTS(componenttol=0.9999)

  So if by “similar component applications” you mean exact local graph-component
  exponentials analogous to Pauli rotations, yes — that is now the core
  backend.

✻ Baked for 21s

❯ no there was a speedup you suggested, where we checked if components could be 
  batched to reduce cost                                                        

● We did not implement batching by component type/signature.

  We only investigated whether batching might be possible.

  What we checked

  We computed component signatures to see whether many components share the same
  structure.

  The script was:

  /home/dhruva/.claude/jobs/513a65b3/tmp/component_signature_stats.jl

  It measured:

  - number of components,
  - component size distribution,
  - “degree signatures”,
  - “spectral signatures” of component matrices,
  - multiplicities of repeated signatures.

  What we found

  There is repetition, but not enough to make batching an obvious easy win.

  Example N=32, M=12:

  ┌───────┬────────────┬──────────────────────────┬─────────────────────────┐
  │ Orbit │ Components │     Unique spectral      │     Unique degree       │
  │       │            │        signatures        │       signatures        │
  ├───────┼────────────┼──────────────────────────┼─────────────────────────┤
  │     1 │        782 │                      121 │                       8 │
  ├───────┼────────────┼──────────────────────────┼─────────────────────────┤
  │     2 │        670 │                      129 │                       9 │
  ├───────┼────────────┼──────────────────────────┼─────────────────────────┤
  │     3 │        697 │                       44 │                       4 │
  └───────┴────────────┴──────────────────────────┴─────────────────────────┘

  Some high-multiplicity signatures existed:

  orbit 3 top spectral multiplicity: 398
  orbit 1 top spectral multiplicities: 179, 112, 109, 92, ...
  orbit 2 top spectral multiplicities: 152, 90, 86, 71, ...

  So batching may be possible, especially for common component types, but we did
  not implement it.

  Why we did not proceed then

  At the time, we first implemented:

  1. per-component exact exponentials,
  2. component plan caching,
  3. component-weight truncation,
  4. freeze/drop/tail approximations.

  Batching was left as a future optimization because it requires a reliable
  canonical component signature and a way to reuse:

  exp(dt * A_signature)

  across components with different node labels but identical matrix structure.

  What batching would mean

  For each component, instead of treating it as unique:

  component C_i -> build A_i -> exp(A_i)

  we would compute a signature/canonical form:

  signature(C_i)

  Then group components:

  all components with same signature share one matrix exponential

  For each component instance, we would only gather/scatter coefficients in its
  local ordering.

  Status

  Not implemented.

  It remains a plausible optimization, and our signature measurements suggest it
  may help — but it is not in the current code.

✻ Worked for 29s

❯ I suppose right now our main issue is that our method is not accurate enough  
  to beat RK. You hypothesized that this is because this is an order 4 method,  
  whereas ours is an order 2 method. How do we study where we can try to        
  improve?                                                                      

● We need to separate the error sources. Right now several things are mixed:

  1. Orbit Strang splitting error: second-order product formula over orbit
  Hamiltonians.
  2. Component truncation/freezing error: we freeze low-weight components.
  3. Post-step trim error: normal coefficient truncation.
  4. RK4 reference error: RK4 is not exact except in small dense tests.
  5. Observable-only error: we mostly compare a small observable set.

  To know where to improve, we need an error budget.

  The core question

  If TrotterTS is inaccurate because of order-2 Strang, then:

  error \sim C_{\mathrm{Strang}} dt^2

  and decreasing dt should reduce error quadratically until another error source
  dominates.

  If it is inaccurate because of component freezing, then decreasing dt will not
  help much; increasing componenttol should help.

  If it is inaccurate because of trim, increasing M should help.

  So we need sweeps over:

  dt
  componenttol
  M

  and compare against exact on small systems.

  ---
  Experiment 1: exact small-system error decomposition

  Use small systems where dense exact is possible:

  N = 4, 6, 8

  Run:

  A. Full component exact, no freeze

  componenttol = 1.0

  This isolates:

  orbit Strang splitting + trim

  B. High componenttol

  componenttol = 0.9999, 0.99999

  This adds freeze error.

  C. No/large trim

  Set M large enough that trim does nothing.

  This isolates pure orbit-Strang error.

  For each:

  dt = 0.2, 0.1, 0.05, 0.025

  Measure exact dense error.

  Expected outcomes:

  - If componenttol=1, large M shows slope ~2 in log-log error vs dt, then
  Strang error dominates.
  - If error plateaus as dt shrinks, then component/trim/backend error
  dominates.
  - If increasing M changes error, trim matters.

  ---
  Experiment 2: compare first-order vs second-order product formula

  TrotterTS(order=1) vs TrotterTS(order=2).

  If second-order is working, error slopes should be:

  order=1: global O(dt)
  order=2: global O(dt^2)

  This validates the orbit product formula implementation itself.

  If order=2 does not show slope ~2, something else is dominating.

  ---
  Experiment 3: freeze-error isolation at fixed dt

  For fixed small dt, compare:

  componenttol = 1.0
  componenttol = 0.999999
  componenttol = 0.9999
  componenttol = 0.999
  componenttol = 0.99

  against exact.

  If componenttol=1.0 is much better than 0.9999, freezing matters.

  If componenttol=1.0 is only slightly better, freezing is not the dominant
  error.

  From earlier sweeps, increasing tolerance above 0.9999 barely improved error,
  suggesting Strang splitting dominates once tol is high enough.

  But we need exact small-system confirmation.

  ---
  Experiment 4: trim-error isolation

  Run exact small tests with:

  M = high enough to avoid trim
  M = normal
  M = aggressive

  If high M and normal M match, trim is not the dominant source.

  ---
  Experiment 5: higher-order product formula over orbit flows

  If Experiment 1 confirms second-order splitting dominates, then the natural
  improvement is a higher-order Suzuki composition over orbit flows.

  For example, fourth-order Suzuki from second-order Strang:

  S_4(dt) =
  S_2(p dt)^2
  S_2((1 - 4p) dt)
  S_2(p dt)^2

  with:

  p = \frac{1}{4 - 4^{1/3}}.

  This uses five second-order steps with scaled timesteps. It costs ~5× more per
  timestep but has global fourth-order accuracy.

  Then compare:

  RK4 dt=0.1
  TrotterTS order-2 dt=0.05 or 0.025
  TrotterTS order-4 dt=0.1

  If order-4 TrotterTS gets RK-level accuracy at fewer steps, it may recover
  speed.

  Caveat: Suzuki fourth order has a negative timestep substep, so component
  freeze/truncation may behave weirdly. But mathematically it is the standard
  path.

  ---
  What I would do next

  Run a small exact diagnostic table:

  N=6 or 8
  H=MFIM
  O=Xtot
  observables = Xtot, Ztot, ZZtot, maybe full matrix norm
  M large
  componenttol = 1.0
  order = 1,2
  dt = 0.2, 0.1, 0.05, 0.025

  Report:

  error
  error ratio when dt halves
  estimated order
  runtime

  This will tell us whether the current method is truly second-order and whether
  splitting error is the main bottleneck.

  Then:

  - if order ~2: implement/test fourth-order Suzuki orbit composition;
  - if error plateaus: investigate component/trim error;
  - if order is wrong: bug in component flow or Strang composition.

  My current hypothesis

  Based on the tolerance sweep:

  componenttol beyond 0.9999 did not improve error much
  dt=0.1 -> 0.05 did improve error a lot

  I expect:

  dominant error = second-order orbit Strang splitting

  So the next likely improvement is a fourth-order orbit product formula, not
  component tolerance tuning.

✻ Sautéed for 35s

❯ Do we need trim, if we are already doing orbit trimming?                      

● Maybe not at every step, and it may be counterproductive to apply both
  aggressively.

  Right now TrotterTS has two truncations:

  1. Orbit-component truncation/freezing inside each orbit flow:

  keep high-weight components exactly
  freeze low-weight components
  2. Normal coefficient trim after each evolution step:

  O = truncation(O)

  2. where in examples:

  truncation(o) = trim(o, 2^M)

  These are different truncation mechanisms.

  What component truncation does

  Component truncation acts at the level of orbit-flow invariant subspaces:

  component C under L_a

  It decides which components to evolve exactly and which to freeze.

  It does not necessarily reduce the number of Pauli strings, because frozen
  components are still carried forward. It reduces expensive work, not
  necessarily final support.

  For freeze mode:

  low-weight components remain in O

  So the operator can still have many terms unless later trim removes them.

  What trim does

  trim(o, M) keeps the largest individual coefficients.

  This controls memory/representation size directly.

  It is independent of orbit components.

  If we remove trim entirely

  Then TrotterTS may keep many strings over time, especially frozen ones, and
  the active basis may grow. Component caching might also grow.

  So no trim could be accurate but slow/memory-heavy.

  But applying trim after every step may be too much

  Since component truncation already approximates the orbit flow, post-step trim
  can introduce an additional error source.

  To isolate error, we should run:

  componenttol = 1.0, no trim
  componenttol = 1.0, trim
  componenttol < 1.0, no trim
  componenttol < 1.0, trim

  on small exact systems.

  This tells us whether normal trim matters for accuracy and support growth.

  Practical possibilities

  Option A: Use only trim

  Set:

  componenttol = 1.0

  so every component is evolved exactly, then use trim.

  This tests pure orbit Trotter + coefficient truncation.

  Option B: Use only component truncation

  Set:

  truncation = identity

  and componenttol < 1.

  This tests whether component truncation alone controls growth.

  Option C: Use both

  Current behavior.

  Option D: Trim less frequently

  Maybe trim only after a full save interval, or every few orbit flows, not
  after every substep. Currently in TrotterTS, truncation is applied after the
  full Trotter step in _evolve, and also inside _orbit_flow only as final
  return:

  return truncation(out)

  Actually, that means each orbit flow applies truncation too. For second-order
  with m=3, that is 5 trims per step plus the outer trim. This may be more than
  needed.

  We should check that. In _orbit_flow, it ends with:

  return truncation(out)

  and _evolve then does:

  O = truncation(O)

  So yes, TrotterTS trims after every orbit flow and again after the full step.

  That could be a significant extra error and cost.

  Immediate recommendation

  For TrotterTS, component truncation is already integrated at the orbit-flow
  level, so I would test moving normal truncation to only after the full Trotter
  step, not after each orbit flow.

  That means changing _orbit_flow to:

  return out

  and let _evolve do:

  O = truncation(O)

  after the full Strang step.

  This may improve accuracy and maybe speed.

  But before changing production behavior, test in a scratch comparison:

  orbit-flow trim on/off
  final-step trim on/off

  Short answer

  No, we may not need ordinary trim inside each orbit flow if component
  truncation is already selecting work. We probably still need some final trim
  to control representation size.

  The best candidate is:

  component truncation inside orbit flow
  normal trim only after full Trotter step

  not both at every orbit subflow.

✻ Crunched for 29s

❯ Ok fine, make the full diagnostic table from the previous comment             

● Write(~/.claude/jobs/513a65b3/tmp/diagnose_trotterts_errors.jl)
  ⎿  Wrote 112 lines to
     ../../../.claude/jobs/513a65b3/tmp/diagnose_trotterts_errors.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +132 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)
● Ran the full diagnostic table to isolate error sources.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/diagnose_trotterts_errors.jl

  Small exact system:

  N=6, M=10, T=1.0

  Observables plus full matrix error were compared against dense exact.

  ---
  1. Order diagnostic, exact components, no trim

  First-order orbit product

  ┌───────┬───────────────┬──────────────────┐
  │  dt   │ obs rel error │ matrix rel error │
  ├───────┼───────────────┼──────────────────┤
  │   0.2 │       9.29e-2 │          1.35e-1 │
  ├───────┼───────────────┼──────────────────┤
  │   0.1 │       4.51e-2 │          6.71e-2 │
  ├───────┼───────────────┼──────────────────┤
  │  0.05 │       2.23e-2 │          3.37e-2 │
  ├───────┼───────────────┼──────────────────┤
  │ 0.025 │       1.11e-2 │          1.69e-2 │
  └───────┴───────────────┴──────────────────┘

  Error ratios when halving dt:

  ~2.06, 2.02, 2.01

  So order 1 behaves correctly.

  Second-order orbit product

  ┌───────┬───────────────┬──────────────────┐
  │  dt   │ obs rel error │ matrix rel error │
  ├───────┼───────────────┼──────────────────┤
  │   0.2 │       2.07e-2 │          3.08e-2 │
  ├───────┼───────────────┼──────────────────┤
  │   0.1 │       5.30e-3 │          7.64e-3 │
  ├───────┼───────────────┼──────────────────┤
  │  0.05 │       1.35e-3 │          1.91e-3 │
  ├───────┼───────────────┼──────────────────┤
  │ 0.025 │       3.42e-4 │          4.76e-4 │
  └───────┴───────────────┴──────────────────┘

  Error ratios when halving dt:

  ~3.90, 3.92, 3.95

  So TrotterTS(order=2, componenttol=1) is genuinely second-order.

  This strongly confirms the dominant high-accuracy error is orbit Strang
  splitting, not implementation bugs.

  ---
  2. Trim isolation

  For componenttol=1, order 2:

  ┌───────┬─────────────────┬──────────────┬──────────────────┐
  │  dt   │ no trim obs rel │ trim obs rel │ RK4 trim obs rel │
  ├───────┼─────────────────┼──────────────┼──────────────────┤
  │   0.1 │        5.304e-3 │     5.304e-3 │         4.435e-4 │
  ├───────┼─────────────────┼──────────────┼──────────────────┤
  │  0.05 │        1.352e-3 │     1.352e-3 │         2.778e-5 │
  ├───────┼─────────────────┼──────────────┼──────────────────┤
  │ 0.025 │        3.420e-4 │     3.420e-4 │         1.736e-6 │
  └───────┴─────────────────┴──────────────┴──────────────────┘

  Normal trim has no effect in this small case because support stays below cap:

  len = 699

  So trim is not the source of the error here.

  ---
  3. Component tolerance isolation

  dt = 0.1

  ┌──────────────┬──────────┬────────────┬─────┐
  │ componenttol │ obs rel  │ matrix rel │ len │
  ├──────────────┼──────────┼────────────┼─────┤
  │         0.99 │ 2.430e-2 │   1.565e-1 │  95 │
  ├──────────────┼──────────┼────────────┼─────┤
  │        0.999 │ 6.592e-3 │   5.231e-2 │ 235 │
  ├──────────────┼──────────┼────────────┼─────┤
  │       0.9999 │ 5.356e-3 │   1.891e-2 │ 427 │
  ├──────────────┼──────────┼────────────┼─────┤
  │      0.99999 │ 5.325e-3 │   9.498e-3 │ 579 │
  ├──────────────┼──────────┼────────────┼─────┤
  │          1.0 │ 5.304e-3 │   7.641e-3 │ 699 │
  └──────────────┴──────────┴────────────┴─────┘

  dt = 0.05

  ┌──────────────┬──────────┬────────────┬─────┐
  │ componenttol │ obs rel  │ matrix rel │ len │
  ├──────────────┼──────────┼────────────┼─────┤
  │         0.99 │ 2.262e-2 │   1.467e-1 │  96 │
  ├──────────────┼──────────┼────────────┼─────┤
  │        0.999 │ 2.452e-3 │   5.122e-2 │ 235 │
  ├──────────────┼──────────┼────────────┼─────┤
  │       0.9999 │ 1.510e-3 │   1.707e-2 │ 452 │
  ├──────────────┼──────────┼────────────┼─────┤
  │      0.99999 │ 1.376e-3 │   5.560e-3 │ 603 │
  ├──────────────┼──────────┼────────────┼─────┤
  │          1.0 │ 1.352e-3 │   1.906e-3 │ 699 │
  └──────────────┴──────────┴────────────┴─────┘

  Interpretation:

  - Component freezing affects matrix/global state error significantly.
  - Observable error plateaus near the exact-component Strang error once
  componenttol >= 0.9999.
  - Raising tolerance improves matrix error, but costs more and increases
  length.

  ---
  4. RK4 baseline order

  ┌───────┬─────────────┬────────────────┐
  │  dt   │ RK4 obs rel │ RK4 matrix rel │
  ├───────┼─────────────┼────────────────┤
  │   0.2 │    6.894e-3 │       2.430e-2 │
  ├───────┼─────────────┼────────────────┤
  │   0.1 │    4.435e-4 │       1.575e-3 │
  ├───────┼─────────────┼────────────────┤
  │  0.05 │    2.778e-5 │       9.879e-5 │
  ├───────┼─────────────┼────────────────┤
  │ 0.025 │    1.736e-6 │       6.178e-6 │
  └───────┴─────────────┴────────────────┘

  RK4 is behaving as fourth order.

  ---
  Main conclusion

  The method is implemented correctly as a second-order orbit product formula.

  The limiting accuracy source at high componenttol is:

  second-order Strang splitting over orbit Hamiltonians

  not:

  trim
  component freezing
  implementation error

  For RK-level accuracy, a second-order method needs much smaller dt, which eats
  speedup.

  Where to improve

  If we want RK4-level accuracy and speedup, the next algorithmic improvement
  should be:

  higher-order orbit product formula

  probably fourth-order Suzuki composition over the existing second-order orbit
  flow.

  Component tolerance tuning alone cannot get RK4-level accuracy at dt=0.1.

✻ Brewed for 2m 11s

❯ Then how is non-TS trotter reaching RK4 ?                                     

● Because in the non-TS benchmark, Trotter is probably not matching RK4 by being
  higher order. It wins because each Trotter substep is much cheaper, so you
  can afford either:

  1. enough small Trotter steps, or
  2. the Trotter error constant is small for that
  Hamiltonian/observable/truncation regime.

  But mathematically, second-order Trotter is still second order. It does not
  generally “reach RK4 accuracy” at the same dt by order alone.

  There are several practical reasons non-TS Trotter can look competitive or
  better.

  ---
  1. Non-TS Trotter applies exact single-Pauli rotations

  For non-TS:

  H = \sum_j c_j P_j.

  Second-order Trotter applies:

  e^{\tau L_{P_j}}

  for each single Pauli term.

  Each one is exact and cheap:

  e^{i\theta P/2} Q e^{-i\theta P/2}
  =
  \cos(\theta) Q + \sin(\theta) Q'.

  So a Trotter step is:

  many very cheap exact Pauli rotations

  RK4 is:

  4 full-H commutators + many operator additions

  Even if Trotter needs smaller dt, it can still win.

  ---
  2. TS orbit Trotter subflows are not cheap single-Pauli rotations

  For TS orbit Trotter:

  H_a^{orb} = \sum_r T_r P_a T_r^\dagger.

  The subflow:

  e^{\tau L_{H_a^{orb}}}

  is not a two-branch Pauli rotation. We implement it via graph components and
  dense exponentials/freeze.

  So TS Trotter has a much more expensive “exact subflow” than non-TS Trotter.

  That is the key difference.

  ---
  3. Current TS RK4 is already TS-native

  In the non-TS case, RK4 operates in ordinary Pauli space.

  In the TS case, RK4 already uses:

  commutator(H::OperatorTS, O::OperatorTS)

  which stays compressed in the TS basis.

  So TS Trotter is competing against a much stronger RK4 baseline than non-TS
  Trotter is.

  This is the big reason the “Trotter beats RK4 non-TS, therefore should beat
  RK4 TS” intuition breaks.

  ---
  4. Non-TS Trotter and TS Trotter split different objects

  Non-TS Trotter splits into single Pauli terms:

  H = \sum_j c_j P_j.

  TS orbit Trotter splits into orbit Hamiltonians:

  H = \sum_a H_a^{orb}.

  Those orbit Hamiltonians are larger, internally nontrivial generators.

  If instead TS Trotter splits into translated individual Pauli terms, then it
  recovers single-Pauli exact rotations — but then it has O(N) gates and loses
  the TS advantage.

  So there is a tradeoff:

  single-Pauli split:
      cheap exact subflows
      many translated gates

  orbit split:
      few factors
      expensive subflows

  The non-TS case gets the best of both: cheap exact subflows and no
  translation-induced gate explosion.

  ---
  5. Error constants may be favorable

  For local Hamiltonians and observables, product formula error can be much
  smaller than pessimistic O(dt^2) estimates. Product formula error depends on
  nested commutators.

  So non-TS Trotter may have a small error constant for the examples, especially
  after truncation.

  But the formal order is still lower than RK4.

  ---
  What to check

  To answer definitively, we should run the same exact-order diagnostic for
  non-TS Trotter:

  For small N, compare dense exact errors:

  RK4 dt = 0.2, 0.1, 0.05, 0.025
  non-TS Trotter order=2 same dt

  Expected:

  RK4 error ratio ~16 when dt halves
  Trotter error ratio ~4 when dt halves

  If non-TS Trotter still “matches” RK4 at practical dt, it is due to small
  constants / cheap cost, not formal order.

  ---
  Bottom line

  Non-TS Trotter reaches useful accuracy because its exact subflows are
  extremely cheap single-Pauli rotations.

  Our TS orbit Trotter is trying to reduce factor count by grouping translations
  into orbit Hamiltonians, but then each subflow becomes a nontrivial
  many-Pauli orbit exponential.

  That is why the speed/accuracy tradeoff is harder in the TS case.

✻ Crunched for 50s

❯ Then how is TS-RK4 able to overcome the O(N) scaling?                         

● TS-RK4 does not eliminate the O(N) relative-translation factor. It overcomes
  the worse scaling by keeping everything in the orbit basis and doing only a
  few such passes.

  The key distinction is:

  TS-RK4 still pays O(N) per commutator,
  but it avoids O(N) separate high-level gate updates and avoids expanding the
  state.

  What TS commutator does

  For:

  commutator(Hts, Ots)

  the TS binary kernel computes:

  for representative h in Hts
    for representative q in Ots
      for relative translation δ
        compute [h, shift(q, δ)]
        fold result back to PauliStringTS

  So the cost is roughly:

  length(Hts) * length(Ots) * N

  For MFIM:

  length(Hts) = 3

  so one full-H TS commutator costs:

  3 * length(Ots) * N

  That is still O(N).

  Why this is still good

  RK4 needs only 4 full-H commutators per time step:

  4 * 3 * length(Ots) * N
  = 12 * length(Ots) * N relative-shift passes

  and each pass immediately folds results back into TS representatives.

  So the state stays compressed.

  What expanded TS Trotter does

  Expanded TS Trotter builds:

  3N translated ordinary gates

  Second-order Strang gives:

  6N - 1 gates

  Each gate updates every current ordinary Pauli string:

  for gate in 1:(6N-1)
      for string in Oordinary
          apply Pauli rotation

  Cost:

  (6N) * length(Oordinary)

  But Oordinary is not the same size as Ots. It is a representative ordinary
  workspace and may carry many translation-equivalent strings that TS-RK4 would
  fold together.

  Also, the gate updates are sequential, with truncation and dictionary rebuilds
  after each gate.

  So even if both have an N factor, RK4 has fewer high-level operations and
  better compression.

  What TS orbit Trotter does

  Our orbit Trotter tries to reduce:

  6N ordinary gates

  to:

  2m - 1 = 5 orbit flows

  But each orbit flow must compute or approximate:

  e^{τ L_{H_a^{orb}}}.

  That is harder than one TS commutator.

  If we use component exponentials, we need component discovery, matrices,
  exponentials, gather/scatter.

  If we use Krylov, we need many commutators and orthogonalization.

  So TS orbit Trotter reduces factor count but increases subflow complexity.

  Why RK4 is hard to beat

  TS-RK4 has the simplest possible use of the TS algebra:

  apply L_H to O four times

  No exact orbit exponentials, no graph components, no many per-component matrix
  exponentials.

  It benefits from TS representation but does not ask for finite orbit flows.

  So it is a strong baseline.

  Simple cost picture

  For MFIM, m=3.

  TS-RK4

  4 full-H commutators
  ≈ 4 * 3 = 12 orbit commutator units

  Orbit Trotter with exact orbit flows

  5 orbit exponentials

  This sounds cheaper, but each orbit exponential is not one commutator unit. It
  may cost many component operations or Krylov applications.

  If each orbit exponential costs more than about:

  12 / 5 ≈ 2.4 orbit-commutator units

  then RK4 wins.

  Our current orbit exponential costs more than that unless we aggressively
  freeze components.

  So how does TS-RK4 “overcome” O(N)?

  It doesn’t remove it. It makes the constant and representation manageable:

  1. Only 4 commutator stages.
  2. Each stage uses TS representatives.
  3. Results fold immediately into OperatorTS.
  4. No sequential O(N) gate list.
  5. No expensive finite orbit exponential machinery.

  That is why TS-RK4 is so competitive.

✻ Sautéed for 29s

❯ So why can't we also reduce the constant in the Expanded TS?                  

● We probably can — and that may be the more practical direction.

  Expanded TS Trotter has the right cheap primitive:

  single-Pauli exact rotations

  but the wrong execution style:

  6N separate gates, ordinary workspace, lots of truncation/dictionary rebuilds

  So instead of trying to exponentiate orbit Hamiltonians, we can try to make
  the expanded translated-gate Trotter execute in a TS-aware way with smaller
  constants.

  This would keep the non-TS Trotter advantage: exact cheap Pauli rotations.

  Where the constant is bad now

  Current TS Trotter does:

  Hr = resum(H)
  gates = trotterize(Hr, dt)
  Or = representative(O)
  trotter_step!(Or, gates; truncation=truncation)

  Problems:

  1. Hr has mN terms.
  2. Strang gives 2mN-1 gates.
  3. trotter_step! applies each gate separately.
  4. Each gate rebuilds an operator dictionary.
  5. Default truncation happens after every gate.
  6. The evolving state is not kept in TS representation during the step.
  7. Translation-equivalent results are not folded until save/final return.

  Even if each gate is cheap, the overhead is huge.

  Possible constant reductions

  1. Lazy translated gates

  We already tried this. It avoids materializing resum(H) but does not reduce
  core work. Minor.

  2. Tune truncate_every

  We found this helps:

  truncate_every = 4–16

  often improved expanded Trotter. But it is a tuning knob, not a structural
  fix.

  Still, it is low-risk and practical.

  3. Fuse gates before rebuilding dictionaries

  Currently each gate does:

  input operator -> dictionary -> output operator

  for every gate.

  Instead, for a block of gates, propagate terms through several gates before
  materializing/truncating.

  But branching grows exponentially unless controlled.

  This is basically what truncate_every approximates.

  4. Fold to TS after blocks of translated gates

  Instead of folding only at the end of the whole evolution, we can do:

  apply some number of translated single-Pauli gates
  fold ordinary workspace to OperatorTS
  continue from representative

  This can reduce duplicate translated equivalents.

  But folding too often may change finite-step ordering/projection, and it costs
  canonicalization.

  We tried folding per full step; it did not help. But folding per block/after
  full Strang step with better truncation semantics might still be worth
  measuring.

  5. Apply translated gate orbits as a batched map

  For a Hamiltonian representative P, instead of treating each translated P_r as
  a separate high-level gate, define an “orbit gate sweep”:

  for q in O:
      apply all translated P_r rotations that affect q

  This keeps single-Pauli rotation physics but batches the orbit.

  For each input string, only translations that overlap/anticommute matter. The
  rest commute and do nothing.

  This is analogous to active-shift orbit_edges, but for finite rotations rather
  than commutators.

  This may be the real middle ground:

  not orbit exponential of sum P_r
  not naive expanded gates
  but active translated single-Pauli rotations in a batched TS-aware sweep

  6. Exploit commutation only generically when detected

  We said not to focus on MFIM structure. But a generic detector could identify
  when translated copies commute and apply them as a batch. However, that is
  extra structure beyond TS alone, so maybe not the core.

  The most promising “constant reduction” path

  I think it is:

  ▎ TS-aware batched translated Pauli rotations.

  For a single translated Pauli gate, ordinary non-TS update is:

  Q \mapsto
  \begin{cases}
  Q, & [P_r,Q]=0 \\
  \cos(\theta) Q + \sin(\theta) Q', & [P_r,Q]\neq0.
  \end{cases}

  For an orbit of translated gates, current expanded Trotter does:

  for r in translations:
      apply rotation P_r

  A batched version would work directly on TS representatives and active
  relative shifts:

  for each q in Ots:
      find active relative shifts r where P_r anticommutes with q
      apply those rotations in the correct Trotter order
      fold outputs to TS representatives

  This preserves the cheap single-Pauli rotation rule while avoiding scanning
  all N translations and avoiding some ordinary workspace duplication.

  Caveat

  Translated gates generally do not commute, so applying all active rotations
  for a string can branch and ordering matters. We must preserve the same
  ordering as expanded Trotter if we want equivalence.

  But for one input TS string q, the number of active shifts may be small. If we
  apply only active shifts, we skip all commuting gates.

  This could reduce the constant substantially while retaining the expanded
  Trotter semantics.

  Why this differs from orbit Liouvillian components

  Orbit Liouvillian components exponentiate the sum:

  e^{τ \sum_r L_{P_r}}

  Expanded Trotter applies product:

  \prod_r e^{τ L_{P_r}}.

  The batched active-gate approach still computes the product formula, but skips
  gates that are no-ops for each branch.

  So it is closer to non-TS Trotter and may preserve its practical speed.

  What to test

  Before implementation, measure for expanded Trotter:

  For each current Pauli string q and each gate orbit P:
      how many translated gates anticommute?

  We already measured commutator active shifts:

  N=64: ~6–7% nonzero
  avg folded edges ~4

  So many gates are no-ops.

  But after branching, active shifts may grow.

  Prototype a single orbit sweep:

  active_orbit_rotation(P, Ots, θ)

  that applies translated single-Pauli rotations only for active shifts, in the
  same order as expanded gates, and compare against explicit expanded orbit
  gates for small systems.

  If equivalent and cheaper, this is likely a better path than orbit
  exponentials.

  Answer

  Yes, we can and probably should reduce constants in expanded TS Trotter.

  The best candidate is not just lazy resum, but skip commuting translated gates
  per string/branch and batch translated Pauli rotations in TS space.

  This preserves the cheap single-Pauli Trotter primitive instead of replacing
  it with expensive orbit exponentials.

✻ Sautéed for 35s

❯ /background                                                                   

※ recap: We’re building a generic TS-native Trotter method that avoids expanded translated-gate overhead while preserving useful accuracy. Next, test a TS-aware batched translated-Pauli rotation that skips commuting/no-op translated gates.

❯ ok. create a branch from main, and try this                                   

  Ran 3 shell commands (ctrl+o to expand)

● Created and switched to a new branch from main:

  ts-active-orbit-trotter

  I’ll proceed there with the active translated single-Pauli orbit-rotation
  prototype.

✻ Sautéed for 1m 35s

❯ Yes go ahead                                                                  

  Ran 1 shell command (ctrl+o to expand)

● The output indicates we're on the new branch from main, so none of the
  previous TrotterTS implementation is present in this branch yet. This is
  expected because we created ts-active-orbit-trotter from main after committing
  the prior work on branch 78.

  Current branch:

  ts-active-orbit-trotter

  Relevant status:

  src/evolution.jl still has only Trotter/RK4/DOPRI5/Exact
  no TrotterTS in this branch

  So before trying active translated single-Pauli orbit sweeps, I need to decide
  whether to:

  1. cherry-pick the previous commits from 78 that added TrotterTS/docs, then
  continue from there, or
  2. implement a narrower prototype directly from main.

  Since this branch was created from main specifically to try a different
  direction, I’ll do option 2: implement only the active
  expanded-gate/orbit-sweep prototype needed for this branch, without bringing
  in the whole previous TrotterTS component machinery.

✻ Cogitated for 44s

❯ yes                                                                           

● Write(~/.claude/jobs/513a65b3/tmp/active_orbit_sweep_prototype.jl)
  ⎿  Wrote 97 lines to
     ../../../.claude/jobs/513a65b3/tmp/active_orbit_sweep_prototype.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 
       4 coord(site,Ls,k) = begin
       5     x=site-1
       6     for j in 1:k-1; x=div(x,Ls[j]); end
       7     mod(x,Ls[k])+1
       8 end
       9 function active_shift(site_g, site_p, Ls, Ps)
      10     valid=true
     … +91 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● I prototyped the active translated single-Pauli orbit sweep from main.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/active_orbit_sweep_prototype.jl

  What it does

  For one Hamiltonian representative G and one ordinary representative workspace
  O, it applies the translated gate orbit:

  for all translations r:
      exp(i θ G_r / 2) O exp(-i θ G_r / 2)

  but it skips gates that commute with the current branch.

  Algorithm per input Pauli string:

  1. Keep a stack of branches:

  (current_pauli, coefficient, last_applied_shift_index)
  2. Find active translated gates after last_applied_shift_index that
  overlap/anticommute.
  3. Apply the earliest active gate.
  4. Branch if noncommuting.
  5. Repeat until no active translated gates remain.
  6. Accumulate result.

  This preserves the same gate ordering as explicit expanded Trotter for that
  orbit sweep, but avoids applying commuting/no-op translated gates.

  Correctness test

  Compared active sweep against explicit expanded translated gates for:

  N = 6
  G = ZZ1111
  O = X1 + Y1Z3 + Z2X4
  theta = 0.07

  Result:

  lenA = 12
  lenE = 12
  err = 0.0

  So the active sweep exactly matched explicit translated single-Pauli gates on
  this test.

  Interpretation

  This is promising because it preserves the cheap exact single-Pauli rotation 
  primitive from non-TS Trotter, rather than replacing it with an orbit
  exponential.

  It is not yet integrated into the package and currently works on ordinary
  Operator workspaces, but it demonstrates the core idea:

  expanded orbit sweep semantics
  but skip commuting translated gates per branch

  Next steps would be:

  1. test more random/local cases,
  2. measure active gate counts vs full N,
  3. integrate as an orbit-sweep helper,
  4. build a second-order TS Trotter that uses these active orbit sweeps instead
  of trotterize(resum(H)).

✻ Brewed for 1m 40s

❯ continue, why do you stop?                                                    

  Read 3 files (ctrl+o to expand)

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line, removed 1 line
      24  export Circuits
      25  export pauli_weight, support, pauli_rotation
      26  export TrotterGate, trotterize, trotter_step!
      27 -export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
         -r, RK4, DOPRI5, Exact                                          
      27 +export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
         +r, TrotterTS, RK4, DOPRI5, Exact                               
      28  
      29  using Random
      30  using LinearAlgebra

● Update(src/evolution.jl)  ⎿  Added 16 lines
      26  Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order)
          , gates)
      27  
      28  """
      29 +    TrotterTS(; order=2, gates=nothing, truncate_every=1)      
      30 +                                                               
      31 +Translation-symmetric product-formula integrator that applies t
         +ranslated single-Pauli                                         
      32 +rotations orbit-by-orbit, skipping translated gates that commut
         +e with each branch.                                            
      33 +This keeps the exact Pauli-rotation primitive of ordinary Trott
         +er while avoiding                                              
      34 +work on no-op translated gates.                                
      35 +"""                                                            
      36 +struct TrotterTS <: AbstractEvolutionMethod                    
      37 +    order::Int                                                 
      38 +    gates::Any                                                 
      39 +    truncate_every::Int                                        
      40 +end                                                            
      41 +TrotterTS(; order::Integer=2, gates=nothing, truncate_every::In
         +teger=1) =                                                     
      42 +    TrotterTS(Int(order), gates, Int(truncate_every))          
      43 +                                                               
      44 +"""                                                            
      45      RK4()
      46  
      47  Classical fixed-step 4th-order Runge–Kutta. Takes one internal 
          step per

● I’ll implement TrotterTS on this clean branch as an active translated
  single-Pauli orbit-sweep method (not orbit exponentials). It will preserve
  expanded Trotter semantics orbit-by-orbit while skipping commuting translated
  gates per branch.

● Update(src/evolution.jl)
  ⎿  Added 1 line, removed 1 line
      127  
      128  # Keyword arguments
      129  - `method::AbstractEvolutionMethod = RK4()`. One of [`Trotter`
           ](@ref),
      130 -  [`RK4`](@ref), [`DOPRI5`](@ref), [`Exact`](@ref).           
      130 +  [`TrotterTS`](@ref), [`RK4`](@ref), [`DOPRI5`](@ref), [`Exac
          +t`](@ref).                                                    
      131  - `truncation`: function `O -> O` applied after every internal
            step. Default
      132    `identity`.
      133  - `dissipation`: function `(O, dt) -> O` applied after every i
           nternal step. The

● Update(src/trotter.jl)
  ⎿  Error editing file

  Read 2 files (ctrl+o to expand)

● Update(src/trotter.jl)
  ⎿  Added 80 lines
      107      end
      108      return O
      109  end
      110 +                                                              
      111 +function _active_shift_tuple(site_g::Integer, site_p::Integer,
          + Ls::Tuple, Ps::Tuple)                                        
      112 +    valid = true                                              
      113 +    shifts = ntuple(length(Ls)) do k                          
      114 +        cg = _coord(site_g, Ls, k)                            
      115 +        cp = _coord(site_p, Ls, k)                            
      116 +        if Ps[k]                                              
      117 +            d = mod(cp - cg, Ls[k])                           
      118 +            iszero(d) ? Ls[k] : d                             
      119 +        else                                                  
      120 +            valid &= cg == cp                                 
      121 +            1                                                 
      122 +        end                                                   
      123 +    end                                                       
      124 +    return valid, shifts                                      
      125 +end                                                           
      126 +                                                              
      127 +function _active_shift_indices(G::PauliString, P::PauliString,
          + Ls::Tuple, Ps::Tuple, shift_index, startidx::Integer)        
      128 +    uG = G.v | G.w                                            
      129 +    uP = P.v | P.w                                            
      130 +    seen = falses(length(shift_index))                        
      131 +    inds = Int[]                                              
      132 +    xG = uG                                                   
      133 +    while !iszero(xG)                                         
      134 +        site_g = trailing_zeros(xG) + 1                       
      135 +        xG &= xG - one(xG)                                    
      136 +        xP = uP                                               
      137 +        while !iszero(xP)                                     
      138 +            site_p = trailing_zeros(xP) + 1                   
      139 +            xP &= xP - one(xP)                                
      140 +            valid, shifts = _active_shift_tuple(site_g, site_p
          +, Ls, Ps)                                                     
      141 +            valid || continue                                 
      142 +            idx = shift_index[shifts]                         
      143 +            idx > startidx || continue                        
      144 +            seen[idx] && continue                             
      145 +            seen[idx] = true                                  
      146 +            push!(inds, idx)                                  
      147 +        end                                                   
      148 +    end                                                       
      149 +    sort!(inds)                                               
      150 +    return inds                                               
      151 +end                                                           
      152 +                                                              
      153 +function _apply_pauli_rotation_branch(G::PauliString, theta::R
          +eal, P::PauliString, c)                                       
      154 +    stheta, ctheta = sincos(theta)                            
      155 +    C, k = commutator(G, P)                                   
      156 +    if k == 0                                                 
      157 +        return ((P, c),)                                      
      158 +    end                                                       
      159 +    phase = (1.0im)^ycount(G)                                 
      160 +    return ((P, c * ctheta), (C, c * (1im * stheta / 2) * phas
          +e * k))                                                       
      161 +end                                                           
      162 +                                                              
      163 +function active_orbit_sweep(O::Operator, G::PauliString, theta
          +::Real, Ls::Tuple, Ps::Tuple; truncation::Function=identity)  
      164 +    shifts = collect(all_shifts(Ls, Ps))                      
      165 +    shift_index = Dict(s => i for (i, s) in enumerate(shifts))
      166 +    shifted_generators = Dict{Int,typeof(G)}()                
      167 +    d = emptydict(O)                                          
      168 +                                                              
      169 +    for (P0, c0) in zip(O.strings, O.coeffs)                  
      170 +        stack = Tuple{typeof(P0),typeof(c0),Int}[(P0, c0, 0)] 
      171 +        while !isempty(stack)                                 
      172 +            P, c, pos = pop!(stack)                           
      173 +            inds = _active_shift_indices(G, P, Ls, Ps, shift_i
          +ndex, pos)                                                    
      174 +            if isempty(inds)                                  
      175 +                setwith!(+, d, P, c)                          
      176 +                continue                                      
      177 +            end                                               
      178 +            idx = first(inds)                                 
      179 +            Gs = get!(shifted_generators, idx) do             
      180 +                shift(G, Ls, Ps, shifts[idx])                 
      181 +            end                                               
      182 +            for (P2, c2) in _apply_pauli_rotation_branch(Gs, t
          +heta, P, c)                                                   
      183 +                iszero(c2) || push!(stack, (P2, c2, idx))     
      184 +            end                                               
      185 +        end                                                   
      186 +    end                                                       
      187 +    out = Operator{keytype(d),valtype(d)}(collect(keys(d)), co
          +llect(values(d)))                                             
      188 +    return truncation(out)                                    
      189 +end                                                           

● Update(src/trotter.jl)
  ⎿  Added 8 lines
      108      return O
      109  end
      110  
      111 +function _coord(site::Integer, Ls::Tuple, dim::Integer)       
      112 +    x = site - 1                                              
      113 +    for k in 1:(dim - 1)                                      
      114 +        x = div(x, Ls[k])                                     
      115 +    end                                                       
      116 +    return mod(x, Ls[dim]) + 1                                
      117 +end                                                           
      118 +                                                              
      119  function _active_shift_tuple(site_g::Integer, site_p::Integer,
            Ls::Tuple, Ps::Tuple)
      120      valid = true
      121      shifts = ntuple(length(Ls)) do k

● Update(src/PauliStrings.jl)
  ⎿  Added 1 line, removed 1 line
      23  export k_local_basis, z_basis, x_basis, y_basis, complete_basis
      24  export Circuits
      25  export pauli_weight, support, pauli_rotation
      26 -export TrotterGate, trotterize, trotter_step!                  
      26 +export TrotterGate, trotterize, trotter_step!, active_orbit_swe
         +ep                                                             
      27  export evolve, EvolutionResult, AbstractEvolutionMethod, Trotte
          r, TrotterTS, RK4, DOPRI5, Exact
      28  
      29  using Random

● Update(src/evolution.jl)
  ⎿  Added 59 lines
      295      return EvolutionResult(collect(tspan), history, OperatorTS
           {Ls,Ps}(Or))
      296  end
      297  
      298 +function _ts_orbit_terms(H::Operator{<:PauliStringTS})        
      299 +    terms = typeof(H)[]                                       
      300 +    for j in 1:length(H)                                      
      301 +        push!(terms, typeof(H)([H.strings[j]], [H.coeffs[j]]))
      302 +    end                                                       
      303 +    return terms                                              
      304 +end                                                           
      305 +                                                              
      306 +function _active_ts_orbit_step(O::Operator, Hterms, dt::Real, 
          +method::TrotterTS, hbar::Real, truncation, Ls, Ps)            
      307 +    method.order ∈ (1, 2) || throw(ArgumentError("order must b
          +e 1 or 2, got $(method.order)"))                              
      308 +    isempty(Hterms) && return O                               
      309 +    function sweep!(Ha, scale)                                
      310 +        for (c, p) in zip(get_coeffs(Ha), Ha.strings)         
      311 +            theta = _trotter_theta(c, dt * scale, hbar, true) 
      312 +            O = active_orbit_sweep(O, representative(p), theta
          +, Ls, Ps; truncation=truncation)                              
      313 +        end                                                   
      314 +    end                                                       
      315 +    if method.order == 1 || length(Hterms) == 1               
      316 +        for Ha in Hterms                                      
      317 +            sweep!(Ha, 1)                                     
      318 +        end                                                   
      319 +    else                                                      
      320 +        L = length(Hterms)                                    
      321 +        for j in 1:(L - 1)                                    
      322 +            sweep!(Hterms[j], 1 / 2)                          
      323 +        end                                                   
      324 +        sweep!(Hterms[L], 1)                                  
      325 +        for j in (L - 1):-1:1                                 
      326 +            sweep!(Hterms[j], 1 / 2)                          
      327 +        end                                                   
      328 +    end                                                       
      329 +    return O                                                  
      330 +end                                                           
      331 +                                                              
      332 +function _evolve(method::TrotterTS, H::Operator{<:PauliStringT
          +S}, O::Operator{<:PauliStringTS}, tspan;                      
      333 +                 truncation, dissipation, fout, hbar)         
      334 +    qubitsize(H) == qubitsize(O) && periodicflags(H) == period
          +icflags(O) ||                                                 
      335 +        throw(DimensionMismatch("H and O must share the same t
          +ranslation-symmetry lattice"))                                
      336 +    Ls = qubitsize(O)                                         
      337 +    Ps = periodicflags(O)                                     
      338 +    n = length(tspan)                                         
      339 +    history = _alloc_history(fout, O, n)                      
      340 +    Hterms = _ts_orbit_terms(H)                               
      341 +    Or = representative(O)                                    
      342 +    for i in ProgressBar(1:(n - 1))                           
      343 +        dt = tspan[i + 1] - tspan[i]                          
      344 +        Or = _active_ts_orbit_step(Or, Hterms, dt, method, hba
          +r, truncation, Ls, Ps)                                        
      345 +        Or = dissipation(Or, dt)                              
      346 +        Or = truncation(Or)                                   
      347 +        _save!(history, fout, OperatorTS{Ls,Ps}(Or), i + 1)   
      348 +    end                                                       
      349 +    return EvolutionResult(collect(tspan), history, OperatorTS
          +{Ls,Ps}(Or))                                                  
      350 +end                                                           
      351 +                                                              
      352 +function _evolve(::TrotterTS, H::AbstractOperator, O::Abstract
          +Operator, tspan;                                              
      353 +                 truncation, dissipation, fout, hbar)         
      354 +    throw(ArgumentError("TrotterTS evolution via `evolve` is i
          +mplemented for `OperatorTS` only, not for `$(typeof(H))`."))  
      355 +end                                                           
      356 +                                                              
      357  function _evolve(::Trotter, H::AbstractOperator, O::AbstractOp
           erator, tspan;
      358                   truncation, dissipation, fout, hbar)
      359      throw(ArgumentError("Trotter evolution via `evolve` is imp
           lemented for `Operator` and " *

  Ran 1 shell command (ctrl+o to expand)

● Next I’ll validate the active-orbit implementation beyond the existing
  evolution smoke tests: exact agreement with explicit translated gates for
  single orbit sweeps, small-system exact accuracy, and runtime/allocations vs
  RK4 and existing expanded Trotter.

● Write(~/.claude/jobs/513a65b3/tmp/validate_active_orbit_trotter.jl)
  ⎿  Wrote 86 lines to
     ../../../.claude/jobs/513a65b3/tmp/validate_active_orbit_trotter.jl
       1 using PauliStrings
       3 using Printf
       4 
       5 function MFIM(N,h)
         ",2; OperatorTS{(N,)}(H)
       7 end
         ; end
       9 function Ztot(N); O=Operator(N); O += "Z",1; OperatorTS{(N,)}(O)
         ; end
         N,)}(O); end
     … +96 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Continued and implemented TrotterTS on this branch as an active translated 
  single-Pauli orbit sweep prototype.

  What was implemented

  Added:

  active_orbit_sweep(O, G, theta, Ls, Ps; truncation=identity)

  This applies a full translated orbit of single-Pauli rotations, but for each
  branch it skips translated gates that commute.

  It preserves expanded Trotter semantics orbit-by-orbit, using the cheap exact
  single-Pauli rotation rule.

  TrotterTS now uses this in the TS evolution path:

  TrotterTS(; order=2, gates=nothing, truncate_every=1)

  Correctness of a single orbit sweep

  Compared active_orbit_sweep against explicit translated gates.

  N=6, G=X11111:   err=0.0, len_active=7,  len_explicit=7
  N=6, G=ZZ1111:   err=0.0, len_active=12, len_explicit=12
  N=8, G=X1111111: err=0.0, len_active=7,  len_explicit=7
  N=8, G=ZZ111111: err=0.0, len_active=12, len_explicit=12

  So for single-orbit sweeps, active skipping matches explicit expanded
  translated gates exactly.

  Method validation / benchmark

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/validate_active_orbit_trotter.jl

  N=6, M=10, T=1, dt=0.1

  ┌──────────────────┬─────────┬─────────┬─────┬──────────┬───────────┐
  │      Method      │  Time   │  Alloc  │ Len │   Norm   │ Exact rel │
  ├──────────────────┼─────────┼─────────┼─────┼──────────┼───────────┤
  │ RK4              │ 0.0092s │ 0.018GB │ 699 │ 0.999426 │   4.44e-4 │
  ├──────────────────┼─────────┼─────────┼─────┼──────────┼───────────┤
  │ Expanded Trotter │ 0.0177s │ 0.027GB │ 391 │ 0.999579 │   3.91e-3 │
  ├──────────────────┼─────────┼─────────┼─────┼──────────┼───────────┤
  │ Active TrotterTS │  0.128s │ 0.236GB │ 396 │ 0.999704 │   3.91e-3 │
  └──────────────────┴─────────┴─────────┴─────┴──────────┴───────────┘

  Active TrotterTS vs RK4:

  rel=3.78e-3

  N=8, M=10

  ┌──────────────────┬─────────┬─────────┬──────┬──────────┬───────────┐
  │      Method      │  Time   │  Alloc  │ Len  │   Norm   │ Exact rel │
  ├──────────────────┼─────────┼─────────┼──────┼──────────┼───────────┤
  │ RK4              │ 0.0166s │ 0.024GB │ 1024 │ 0.999426 │   4.44e-4 │
  ├──────────────────┼─────────┼─────────┼──────┼──────────┼───────────┤
  │ Expanded Trotter │ 0.0209s │ 0.032GB │  425 │ 0.998555 │   1.21e-3 │
  ├──────────────────┼─────────┼─────────┼──────┼──────────┼───────────┤
  │ Active TrotterTS │  0.136s │ 0.243GB │  422 │ 0.999738 │   5.32e-3 │
  └──────────────────┴─────────┴─────────┴──────┴──────────┴───────────┘

  N=32, M=12

  ┌──────────────────┬────────┬─────────┬──────┬──────────┐
  │      Method      │  Time  │  Alloc  │ Len  │   Norm   │
  ├──────────────────┼────────┼─────────┼──────┼──────────┤
  │ RK4              │ 0.127s │ 0.115GB │ 4096 │ 0.999426 │
  ├──────────────────┼────────┼─────────┼──────┼──────────┤
  │ Expanded Trotter │ 0.302s │ 0.449GB │ 1452 │ 0.998870 │
  ├──────────────────┼────────┼─────────┼──────┼──────────┤
  │ Active TrotterTS │  1.10s │ 1.557GB │ 1451 │  1.00001 │
  └──────────────────┴────────┴─────────┴──────┴──────────┘

  Active TrotterTS vs RK4:

  rel=3.79e-3

  N=64, M=14

  ┌──────────────────┬────────┬─────────┬───────┬──────────┐
  │      Method      │  Time  │  Alloc  │  Len  │   Norm   │
  ├──────────────────┼────────┼─────────┼───────┼──────────┤
  │ RK4              │ 0.837s │ 0.536GB │ 16384 │ 0.999426 │
  ├──────────────────┼────────┼─────────┼───────┼──────────┤
  │ Expanded Trotter │  2.24s │ 3.476GB │  4993 │ 0.999546 │
  ├──────────────────┼────────┼─────────┼───────┼──────────┤
  │ Active TrotterTS │  6.17s │ 8.053GB │  4998 │  1.00000 │
  └──────────────────┴────────┴─────────┴───────┴──────────┘

  Active TrotterTS vs RK4:

  rel=3.65e-3

  Interpretation

  The active-orbit sweep is correct but currently slower than expanded Trotter.

  Why:

  - It avoids no-op gates, but it uses a branch stack per input term.
  - It recomputes active shifts per branch.
  - It allocates heavily (Dict, Vector, branch tuples).
  - Expanded Trotter’s low-level loop is much tighter, even though it does more
  conceptual gate checks.

  So this direction only becomes viable if implemented much more tightly, likely
  with fewer per-branch allocations and perhaps block-level truncation. As
  written, it verifies the idea but is not performant.

✻ Churned for 7m 51s

❯ explain to me what the method actually is                                     

● The method we just tried is not the orbit-Liouvillian/component-exponential
  method from the previous branch.

  This branch’s method is:

  ▎ Keep the ordinary single-Pauli Trotter rotation rule, but apply each 
  ▎ translated Hamiltonian orbit in a TS-aware “active sweep” that skips 
  ▎ translated gates that commute with the current Pauli branch.

  It is closer to expanded Trotter than to orbit-exponential Trotter.

  ---
  1. Starting point: current expanded TS Trotter

  For a TS Hamiltonian term:

  H_a^{TS} = c_a \sum_r T_r P_a T_r^\dagger

  the current implementation expands it into all translated ordinary Pauli
  terms:

  c_a P_{a,r}.

  Then Strang Trotter applies ordinary single-Pauli rotations:

  e^{\tau L_{P_{a,r}}}

  for every translation r.

  So one orbit sweep is basically:

  for r in all translations:
      apply ordinary Pauli rotation generated by P_a shifted by r

  This uses the cheap non-TS primitive:

  Q \mapsto
  \begin{cases}
  Q, & [P,Q]=0,\\
  \cos(\theta)Q + \sin(\theta)Q', & [P,Q]\neq 0.
  \end{cases}

  The issue is that most translated generators commute with a given Pauli string
  and do nothing.

  ---
  2. Active orbit sweep idea

  Instead of blindly applying every translated gate:

  r = 1, 2, ..., N

  we take one input Pauli string Q and ask:

  ▎ Which translated copies of P_a actually anticommute with Q?

  Only those can change Q.

  So for one branch, we compute the active translated gates:

  active shifts = { r : [P_a shifted by r, Q] ≠ 0 }

  Then we apply only the earliest active gate in the product order.

  If that gate branches:

  Q \to \cos(\theta)Q + \sin(\theta)Q',

  then each branch gets its own remaining active gates, because Q' may
  anticommute with different translated gates than Q.

  So the method is recursive/stack-based.

  ---
  3. Algorithm for one orbit sweep

  Given:

  O = sum c_i Q_i
  G = representative Hamiltonian Pauli string
  theta = rotation angle

  For each Pauli string Q_i in O:

  1. Start a stack with:

  (Q_i, c_i, last_applied_shift = 0)
  2. Pop a branch:

  (Q, c, pos)
  3. Find translated copies of G after pos that anticommute with Q.
  4. If none exist:

  accumulate c * Q into output
  5. Otherwise, take the earliest active translated gate G_r.
  6. Apply the exact ordinary Pauli rotation:

  Q \mapsto \cos(\theta)Q + \sin(\theta)Q'
  7. Push resulting branches back onto the stack with:

  last_applied_shift = r
  8. Repeat until stack is empty.

  This produces exactly the same result as applying all translated ordinary
  gates in order, because commuting gates would have done nothing anyway.

  ---
  4. Why this is still Trotter

  It preserves the same finite product formula as expanded Trotter.

  For one orbit:

  \prod_r e^{\tau L_{P_r}}

  we are not replacing this by:

  e^{\tau \sum_r L_{P_r}}.

  We are still applying the product over translated single-Pauli rotations. We
  are just skipping factors that are the identity on the current branch.

  So this method is:

  active expanded-orbit Trotter

  not:

  orbit-exponential Trotter

  ---
  5. How second-order is built

  For the TS Hamiltonian with orbit terms:

  H = \sum_a H_a^{TS},

  second-order Strang over orbits is:

  half sweep orbit 1
  half sweep orbit 2
  ...
  full sweep orbit m
  ...
  half sweep orbit 2
  half sweep orbit 1

  Each “sweep orbit a” is the active translated single-Pauli gate sweep
  described above.

  This matches the structure of expanded Trotter grouped by orbit, but avoids
  many no-op gates within each orbit sweep.

  ---
  6. What it is supposed to improve

  Expanded Trotter does:

  for every translated gate:
      for every string:
          check commutator

  The active method does:

  for every string/branch:
      find only translated gates that can affect it

  So if most translated gates commute, it can reduce the number of actual
  rotation operations.

  For local strings in large systems, only a small fraction of translations
  overlap/anticommute.

  ---
  7. Why it was slow in the prototype

  The prototype had the right semantics but a poor implementation.

  It did a lot of per-branch work:

  - recomputed active shifts repeatedly,
  - allocated stacks and tuples,
  - used dictionaries,
  - sorted active shift lists,
  - cached shifted generators in small dictionaries.

  So even though it skipped conceptual no-op gates, it was slower than the tight
  existing expanded Trotter loop.

  The expanded loop does more checks, but it is simple and cache-friendly.

  ---
  8. Difference from previous component method

  Previous branch:

  orbit Liouvillian component exponential

  applies:

  e^{\tau L_{H_a^{TS}}}
  =
  e^{\tau \sum_r L_{P_r}}.

  This branch:

  active orbit sweep

  applies:

  \prod_r e^{\tau L_{P_r}}.

  So:

  - component method = exact flow of orbit Hamiltonian, then Strang over orbits;
  - active sweep method = product formula over translated gates, but skipping
  no-op gates.

  They are different approximations.

  ---
  9. Why this direction is attractive

  It keeps the biggest advantage of ordinary non-TS Trotter:

  exact single-Pauli rotations are cheap

  It avoids the hard problem of exponentiating an entire orbit Hamiltonian.

  If implemented efficiently, it could reduce constants without changing the
  basic Trotter semantics.

  ---
  10. Current status

  We verified for single orbit sweeps:

  active sweep == explicit expanded translated gates

  with zero error on test cases.

  But the current implementation is slower than expanded Trotter due to
  overhead.

  So this is a correct prototype, not yet a performant implementation.

✻ Cogitated for 35s

❯ ok why can't we just apply the cheap paulistring update directly to the TS    
  representation, without exanding it out? How is RK able to manage this?       

● Because the cheap Pauli rotation update is not translation invariant at the 
  level of one translated gate.

  RK works because the commutator with the full TS Hamiltonian is translation 
  invariant. A single translated Pauli rotation is not.

  That is the key difference.

  ---
  1. What the cheap Pauli update does

  For an ordinary single Pauli generator G, the update is:

  Q \mapsto
  \begin{cases}
  Q, & [G,Q]=0, \\
  \cos(\theta) Q + \sin(\theta) Q', & [G,Q]\neq 0.
  \end{cases}

  This is cheap because G is a single Pauli string.

  But in the TS Hamiltonian, the generator is not just one G. It is all
  translations:

  [G] = \sum_r G_r.

  The expanded Trotter product uses individual gates:

  e^{\theta L_{G_1}}
  e^{\theta L_{G_2}}
  \cdots
  e^{\theta L_{G_N}}.

  Each factor e^{θ L_{G_r}} singles out an origin/site/bond r.

  That single gate does not preserve translation symmetry.

  ---
  2. Why plain OperatorTS cannot represent intermediate states

  OperatorTS stores translation orbits:

  [Q] = sum_x T_x Q T_x†

  It forgets absolute position.

  But after applying one translated gate G_r, the result depends on the relative
  position of Q to that particular r.

  Example:

  apply gate at site 1

  to a TS orbit [Q].

  The component of [Q] whose representative overlaps site 1 changes. The
  translates of Q far away may not change. So the intermediate state is no
  longer translation invariant.

  Expanded Trotter’s intermediate states are generally not TS after each
  individual translated gate.

  Only after applying all translations in a symmetric/full way do we get
  something close to TS again.

  So if we force-fold after every gate, we lose information about which part of
  the orbit has already interacted with which translated gates.

  ---
  3. How RK avoids this

  RK does not apply individual translated gates.

  It uses the full TS Hamiltonian commutator:

  L_H(O) = i[H^{TS}, O].

  This operator is translation equivariant:

  L_H(TOT†) = T L_H(O) T†.

  Therefore, if O is translation symmetric, L_H(O) is also translation
  symmetric.

  So RK can stay inside OperatorTS:

  OperatorTS -> commutator(Hts, Ots) -> OperatorTS

  It never needs to represent non-TS intermediate states caused by a single
  translated gate.

  ---
  4. Trotter has two choices

  Choice A: Apply individual translated gates

  This keeps the cheap single-Pauli update:

  G_r rotation is cheap

  but intermediate states are not TS.

  So you need either:

  - ordinary Operator workspace,
  - or a richer representation that tracks origin/anchor information.

  Plain OperatorTS is insufficient.

  Choice B: Apply whole orbit flow

  This preserves TS:

  e^{\tau L_{[G]}}

  maps OperatorTS -> OperatorTS.

  But the generator [G] is a sum of translated Paulis, not a single Pauli, so
  the cheap two-branch Pauli rotation rule no longer applies directly.

  Then you need:

  - orbit graph component exponentials,
  - Krylov,
  - Taylor,
  - etc.

  That is what we explored.

  ---
  5. Why active sweep is hard

  The active sweep tried to keep choice A but skip no-op translated gates.

  However, because intermediate branches are not TS, it used an ordinary
  representative workspace and branch stack.

  A truly TS-native version of active sweep would need to remember relative
  origin/phase information that OperatorTS discards.

  That would be an anchored TS representation, e.g.

  (canonical orbit representative, anchor position / sweep cursor) ->
  coefficient

  This would keep enough data to apply individual translated gates without fully
  expanding to ordinary Pauli strings.

  But plain OperatorTS alone cannot do it.

  ---
  6. Summary

  RK can update OperatorTS directly because:

  full TS commutator preserves translation symmetry

  Cheap Pauli rotations cannot be applied directly to OperatorTS because:

  single translated gate breaks translation symmetry

  To use cheap Pauli rotations, we need a representation that can hold non-TS
  intermediate states, either:

  ordinary Operator

  or

  anchored/relative TS representation

  Plain OperatorTS is too compressed; it has quotiented away the origin data
  needed during the ordered product.

✻ Baked for 33s

❯ I see so orbit splitting allows us to apply the trotterization at the orbit   
  level rather than at the full operator level, but the orbit update can be     
  does in a translation symmetric way?                                          

● Yes — exactly.

  Orbit splitting changes the Trotter factors from individual translated Pauli 
  gates to translation-invariant orbit Hamiltonian flows.

  That is the key conceptual move.

  Expanded translated-gate Trotter

  Current expanded TS Trotter splits the Hamiltonian like:

  H = \sum_{a,r} c_a P_{a,r}

  and applies factors:

  e^{\tau L_{P_{a,r}}}.

  Each P_{a,r} is one translated Pauli string.

  These factors are cheap, but not translation symmetric:

  OperatorTS -> non-TS intermediate

  So plain OperatorTS cannot represent the intermediate evolution exactly.

  Orbit-level Trotter

  Orbit splitting instead writes:

  H = \sum_a H_a^{orb},
  \qquad
  H_a^{orb} = c_a \sum_r P_{a,r}.

  Then Trotter applies factors:

  e^{\tau L_{H_a^{orb}}}.

  Each H_a^{orb} is translation symmetric.

  Therefore:

  OperatorTS -> OperatorTS

  for each orbit flow.

  So yes: the orbit update can be done in a translation-symmetric way.

  The tradeoff

  The orbit flow is TS-preserving, but it is no longer a cheap single-Pauli
  rotation.

  For a single Pauli gate:

  e^{\tau L_P}

  the update of one Pauli string is at most two branches.

  For an orbit Hamiltonian:

  e^{\tau L_{\sum_r P_r}},

  the update couples a Pauli orbit to a graph of other Pauli orbits.

  So the TS-native orbit update is:

  build/apply the orbit Liouvillian graph
  exponentiate the connected component or approximate it

  That is exactly what the component method was doing.

  Why this is mathematically clean

  Because each orbit flow is generated by a TS operator:

  L_a(O) = i[H_a^{orb}, O].

  If O is TS, then L_a(O) is TS, and so:

  e^{\tau L_a} O

  is TS.

  Thus the full orbit-level Strang step:

  e^{dt/2 L_1}
  e^{dt/2 L_2}
  ...
  e^{dt L_m}
  ...
  e^{dt/2 L_1}

  stays in OperatorTS at every substep.

  Summary

  Yes:

  orbit splitting = Trotterization at the translation-orbit level
  orbit update = TS-preserving finite flow

  The hard part is that the orbit update is not the ordinary cheap Pauli
  rotation; it is a finite-time evolution on the Pauli-orbit graph.

✻ Crunched for 25s

❯ Ok, can you explain to me how what orbit components are?                      

● An orbit component is the TS analogue of the little two-state subspace that
  appears in an ordinary Pauli rotation.

  Let me build it up from the non-TS case.

  ---
  1. Ordinary Pauli rotation: tiny component

  Suppose we rotate a Pauli string Q by a single Pauli generator P.

  The relevant Liouvillian is:

  L_P(Q) = i[P, Q].

  There are only two cases.

  If P commutes with Q:

  L_P(Q) = 0

  so Q is alone in its component.

  If P anticommutes with Q, then:

  L_P(Q) \propto P Q.

  And applying L_P again sends you back to Q:

  L_P(PQ) \propto Q.

  So the closed subspace is:

  { Q, P Q }

  That is a component of size 2.

  The Pauli rotation is just the exponential of a 2×2 matrix on this component,
  giving the familiar cos/sin formula.

  ---
  2. TS Pauli orbit basis

  In OperatorTS, we do not store a single Pauli string Q. We store its
  translation orbit:

  [Q] = \sum_x T_x Q T_x^\dagger.

  So a basis element is not:

  Q

  but:

  [Q]

  where [Q] means all translations of Q, folded into one canonical
  representative.

  ---
  3. One orbit Hamiltonian

  A TS Hamiltonian representative P defines an orbit Hamiltonian:

  [P] = \sum_r T_r P T_r^\dagger.

  Its Liouvillian is:

  L_{[P]}([Q]) = i[[P], [Q]].

  When you apply this to [Q], you get a linear combination of other TS Pauli
  orbits:

  L_{[P]}([Q]) =
  a_1 [R_1] + a_2 [R_2] + \cdots.

  These [R_i] are found by placing translated copies of P relative to Q, taking
  Pauli commutators, and folding the results back into TS representatives.

  ---
  4. Orbit graph

  This defines a graph:

  node = TS Pauli orbit [Q]
  edge = nonzero coupling under L_[P]

  So if:

  L_{[P]}([Q]) = 2i [R] - 3i [S],

  then the graph has edges:

  [Q] -> [R]
  [Q] -> [S]

  The edge weights are the coefficients.

  ---
  5. Orbit component

  The orbit component of [Q] is the connected/closed set of all TS Pauli orbits
  reachable by repeatedly applying L_[P].

  Start from [Q]:

  C = { [Q] }

  Apply L_[P]:

  [Q] -> [R1], [R2], ...

  Add those.

  Then apply L_[P] to each new one:

  [R1] -> ...
  [R2] -> ...

  Keep going until no new TS Pauli orbits appear.

  That final closed set is the orbit component:

  C([Q]) = { [Q], [R1], [R2], ... }

  It is the smallest subspace such that:

  L_{[P]}(C) \subseteq C.

  ---
  6. Why components matter

  Because on one component, the orbit Liouvillian is just a finite matrix:

  A_C.

  Then the exact orbit flow is:

  e^{\tau L_{[P]}}

  restricted to that component:

  c_C(\tau) = e^{\tau A_C} c_C(0).

  So the component is where the exact TS orbit rotation happens.

  ---
  7. Ordinary Pauli rotation is a special case

  For one ordinary Pauli generator:

  component size = 1 or 2

  So the exponential is tiny and gives cos/sin.

  For one TS orbit generator:

  component size may be 4, 8, 16, 32, ...

  So the exponential is larger, but conceptually the same thing.

  That is why we said:

  ordinary Pauli rotation = exact exponential on a 1- or 2-node component
  TS orbit rotation       = exact exponential on a Pauli-orbit component

  ---
  8. Concrete code picture

  The component is built with repeated calls to:

  orbit_edges(Ha, q)

  where:

  - Ha is one Hamiltonian orbit term,
  - q is one PauliStringTS.

  Pseudo-code:

  component = Set([q])
  frontier = [q]

  while !isempty(frontier)
      q = pop!(frontier)
      edges = orbit_edges(Ha, q)

      for r in edges.strings
          if r not in component
              push!(component, r)
              push!(frontier, r)
          end
      end
  end

  Then we build the component matrix:

  A[i, j] = coefficient of component[i] in L_Ha(component[j])

  and apply:

  exp(dt * A)

  to the coefficient vector on that component.

  ---
  9. Why many components appear

  The active operator O is a sum of many TS Pauli orbits:

  O = \sum_j c_j [Q_j].

  Each [Q_j] belongs to some component under each Hamiltonian orbit Ha.

  So for one orbit flow, the support of O may decompose into many disconnected
  components:

  component 1: [Q1], [R1], ...
  component 2: [Q7], [R9], ...
  component 3: ...

  These components do not talk to each other under this Ha, so they can be
  evolved independently.

  That is why component-wise evolution is mathematically valid.

  ---
  10. The practical issue

  If there are only a few components, this is great.

  If there are hundreds or thousands of components, then even though each is
  small, the overhead of:

  building components
  building matrices
  exponentiating matrices
  gathering/scattering coefficients

  can dominate.

  That is the performance problem we observed.

  ---
  Short version

  An orbit component is:

  ▎ the closed set of translation-symmetric Pauli orbits coupled together by the
  ▎ commutator with one translation-symmetric Hamiltonian orbit.

  It is the TS generalization of the {Q, P Q} pair in an ordinary Pauli
  rotation.

✻ Churned for 39s
❯ If I remember correctly the bulk of it was the fact that there were so many   
  components to update -  which is why we introduced freezing - and not the     
  other steps as much.                                                          

● Yes, that’s right.

  For the component-exponential backend, the main practical problem was not that
  each individual component was huge. Most components were small.

  The problem was:

  there were too many components to update

  Each component update has overhead:

  find/gather active coefficients
  maybe build or look up component plan
  matrix-vector multiply
  scatter back into OperatorTS

  Even if the component is only size 4 or 8, doing that hundreds or thousands of
  times per orbit flow becomes expensive.

  What we measured

  For example, for an evolved N=32, M=12 state:

  orbit 1: ~782 components
  orbit 2: ~670 components
  orbit 3: ~697 components

  So one full Strang step touches something like:

  5 orbit flows × hundreds of components

  That is the bottleneck.

  The component sizes were mostly modest:

  p50 ~ 4–8
  p90 ~ 16–32

  So the bottleneck was component count, not component size.

  Why freezing helped

  The weight distribution was extremely concentrated.

  For N=32:

  ~4–5 components carried 90% of component weight
  ~8–13 components carried 99%
  ~19–24 components carried 99.9%
  ~30–38 components carried 99.99%

  out of hundreds of components.

  So freezing let us avoid expensive component exponentials for most components:

  high-weight components: exact exp update
  low-weight components: carry forward unchanged

  This reduced the number of components we actively updated.

  What freezing did not fix

  Freezing still has to identify components and assign weights, but the costly
  part — dense exponential application and coefficient scatter for every
  component — is reduced.

  And accuracy is limited because frozen components are not evolved.

  The failed tail update

  When we replaced freezing with a cheap first-order update for the tail,
  runtime and allocations grew again because we were once again
  touching/updating all components.

  So yes: the problem was primarily the number of components being actively
  updated. Freeze worked because it skipped most updates.

✻ Brewed for 26s

※ recap: We’re trying to make TS Trotter faster by applying active translated Pauli orbit sweeps instead of expanding every translated gate. Next, continue optimizing the active sweep so it skips no-op gates without heavy per-branch overhead.

❯ The hurdle we faced with the orbit evolution was that the accuracy was not    
  matching rk4. However, our accuracy should atleast match the expanded TS?     

● Yes — that’s a key point.

  At the same product-formula order and timestep, a TS-native Trotter should at
  least be compared to, and ideally match, the expanded TS Trotter, not RK4.

  But there is an important distinction between the two TS-native strategies we
  tried:

  1. Orbit-exponential Trotter does not match expanded TS Trotter exactly

  Orbit-level Trotter applies:

  e^{\tau L_{\sum_r P_r}}

  for each orbit Hamiltonian.

  Expanded TS Trotter applies:

  \prod_r e^{\tau L_{P_r}}.

  These are not the same unless the translated copies commute, or in the limit τ
  -> 0.

  So if by “orbit evolution” you mean the component-exponential method, then no,
  it should not exactly match expanded TS Trotter at finite dt.

  It is a different product formula:

  Strang over orbit Hamiltonians

  rather than:

  Strang over all translated Pauli terms

  Both are second-order approximations to the same exact flow, but with
  different error constants.

  So its accuracy can be worse or better than expanded TS Trotter depending on
  commutators.

  2. Active translated-gate sweep should match expanded TS Trotter

  The active-orbit sweep method does preserve the expanded product formula
  semantics.

  It applies:

  \prod_r e^{\tau L_{P_r}}

  but skips r where the gate commutes with the current branch.

  So this method should match expanded TS Trotter exactly, up to:

  - truncation differences,
  - ordering differences,
  - floating-point accumulation differences.

  We verified this for individual orbit sweeps:

  active_orbit_sweep == explicit translated gates

  with zero error on small tests.

  But the full TrotterTS active implementation may differ from expanded TS
  Trotter if:

  1. The orbit ordering differs from the resum(H) ordering.
  2. Truncation is applied at different points.
  3. The state is folded or represented differently between substeps.
  4. The active sweep applies translations in a different order than
  trotterize(resum(H)).

  What we should test

  We need to compare:

  expanded TS Trotter
  active-orbit TrotterTS

  under identical conditions:

  - same Hamiltonian,
  - same dt,
  - same initial state,
  - no truncation or identical truncation points,
  - same orbit/translation ordering,
  - compare full final Operator or OperatorTS, not just observables.

  For small systems:

  expanded = evolve(Hts, Ots; method=Trotter(), truncation=identity)
  active   = evolve(Hts, Ots; method=TrotterTS(), truncation=identity)
  norm(resum(expanded.final) - resum(active.final))

  If this is not near zero, then active ordering/truncation semantics differ.

  Expected result

  - Component-exponential orbit Trotter: should not necessarily match expanded
  TS.
  - Active translated-gate sweep Trotter: should match expanded TS if
  implemented/order-aligned correctly.

  Since the current branch is active translated-gate sweep, yes: it should at
  least match expanded TS Trotter when truncation/order are aligned.

  That is the next thing to verify.

✻ Churned for 30s

❯ ok let us stop thinking about active translated-gate sweep, and go back to    
  orbit evolution. Switch branches back to 78                                   

  Ran 2 shell commands (ctrl+o to expand)

● Switched back to branch 78.

  I stashed the active translated-gate sweep prototype first, so it is not lost:

  stash@{0}: On ts-active-orbit-trotter: active-orbit-sweep-prototype

  Current branch/status:

  ## 78...origin/78

  Working tree is clean.

✻ Brewed for 1m 9s

❯ First test TrotterTS, Trotter and RK4 for variables. Then, we will explore    
  improving TrotterTS accuracy                                                  

  Searched for 1 pattern (ctrl+o to expand)

● Write(~/.claude/jobs/513a65b3/tmp/compare_methods_variables.jl)
  ⎿  Wrote 97 lines to
     ../../../.claude/jobs/513a65b3/tmp/compare_methods_variables.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       4 
       5 function MFIM(N,h)
       6     H=Operator(N)
       7     H += -h,"X",1
       8     H += -h/2,"Z",2
       9     H += "Z",1,"Z",2
      10     OperatorTS{(N,)}(H)
     … +100 lines (ctrl+o to expand)

● I’ll run a variable sweep comparing RK4, current expanded Trotter, and current
  component-based TrotterTS, including time, allocations, final length, norm,
  and errors vs RK4 / exact where possible.

  Ran 2 shell commands (ctrl+o to expand)

● Here are the variable comparison results for branch 78 current TrotterTS,
  existing expanded Trotter, and RK4.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/compare_methods_variables.jl

  Observables:

  Xtot, Ztot, ZZtot, YZtot

  Small cases also compare to dense exact.

  ---
  N=6, M=10, dt=0.1, T=1.0

  ┌───────────┬──────────┬─────────┬───────────┬──────────┐
  │  Method   │   Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ RK4       │ 0.00939s │ 0.018GB │       699 │ 0.999426 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ Trotter   │  0.0181s │ 0.027GB │       409 │  1.00004 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │  0.0199s │ 0.069GB │       419 │      1.0 │
  └───────────┴──────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4   rel=3.56e-3, max=1.66e-2
  TrotterTS vs RK4   rel=5.35e-3, max=2.61e-2

  RK4       vs exact rel=4.44e-4, max=3.14e-3
  Trotter   vs exact rel=3.71e-3, max=1.83e-2
  TrotterTS vs exact rel=5.34e-3, max=2.61e-2

  ---
  N=8, M=10, dt=0.1, T=1.0

  ┌───────────┬─────────┬─────────┬───────────┬──────────┐
  │  Method   │  Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ RK4       │ 0.0256s │ 0.024GB │      1024 │ 0.999426 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ Trotter   │ 0.0208s │ 0.032GB │       415 │ 0.999979 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.0220s │ 0.073GB │       472 │      1.0 │
  └───────────┴─────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4   rel=9.14e-4, max=6.92e-3
  TrotterTS vs RK4   rel=5.37e-3, max=3.45e-2

  RK4       vs exact rel=4.44e-4, max=4.18e-3
  Trotter   vs exact rel=1.08e-3, max=7.13e-3
  TrotterTS vs exact rel=5.36e-3, max=3.46e-2

  ---
  N=8, M=10, dt=0.05, T=1.0

  ┌───────────┬─────────┬─────────┬───────────┬──────────┐
  │  Method   │  Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ RK4       │ 0.0401s │ 0.057GB │      1024 │ 0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ Trotter   │ 0.0461s │ 0.070GB │       402 │ 0.999904 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.0496s │ 0.153GB │       505 │      1.0 │
  └───────────┴─────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4   rel=2.50e-4, max=1.66e-3
  TrotterTS vs RK4   rel=1.51e-3, max=8.40e-3

  RK4       vs exact rel=2.78e-5, max=2.50e-4
  Trotter   vs exact rel=2.63e-4, max=1.71e-3
  TrotterTS vs exact rel=1.51e-3, max=8.40e-3

  ---
  N=16, M=10, dt=0.1, T=1.0, tol=0.999

  ┌───────────┬──────────┬─────────┬───────────┬──────────┐
  │  Method   │   Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ RK4       │  0.0218s │ 0.027GB │      1024 │ 0.999426 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ Trotter   │  0.0306s │ 0.056GB │       414 │  1.00006 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.00860s │ 0.024GB │       195 │      1.0 │
  └───────────┴──────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4 rel=3.70e-3, max=4.49e-2
  TrotterTS vs RK4 rel=5.92e-3, max=8.19e-2

  ---
  N=16, M=10, dt=0.05, T=1.0, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬──────────┐
  │  Method   │  Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ RK4       │ 0.0495s │ 0.063GB │      1024 │ 0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ Trotter   │ 0.0680s │ 0.119GB │       403 │  1.00003 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.0527s │ 0.186GB │       497 │      1.0 │
  └───────────┴─────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4 rel=9.71e-4, max=1.22e-2
  TrotterTS vs RK4 rel=1.30e-3, max=1.75e-2

  ---
  N=32, M=12, dt=0.1, T=1.0, tol=0.999

  ┌───────────┬──────────┬─────────┬───────────┬──────────┐
  │  Method   │   Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ RK4       │   0.133s │ 0.118GB │      4096 │ 0.999426 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ Trotter   │   0.242s │ 0.426GB │      1403 │  1.00003 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.00877s │ 0.029GB │       239 │      1.0 │
  └───────────┴──────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4 rel=5.66e-3, max=1.54e-1
  TrotterTS vs RK4 rel=6.63e-3, max=1.47e-1

  ---
  N=32, M=12, dt=0.05, T=1.0, tol=0.9999

  ┌───────────┬─────────┬─────────┬───────────┬──────────┐
  │  Method   │  Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ RK4       │  0.295s │ 0.288GB │      4096 │ 0.999982 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ Trotter   │  0.539s │ 0.917GB │      1360 │  1.00001 │
  ├───────────┼─────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.0557s │ 0.206GB │       503 │      1.0 │
  └───────────┴─────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4 rel=1.44e-3, max=3.85e-2
  TrotterTS vs RK4 rel=1.53e-3, max=3.30e-2

  ---
  N=64, M=14, dt=0.1, T=0.5, tol=0.999

  ┌───────────┬──────────┬─────────┬───────────┬──────────┐
  │  Method   │   Time   │  Alloc  │ Final len │   Norm   │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ RK4       │  0.0876s │ 0.077GB │     16384 │ 0.999713 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ Trotter   │   0.687s │ 1.330GB │      4519 │      1.0 │
  ├───────────┼──────────┼─────────┼───────────┼──────────┤
  │ TrotterTS │ 0.00110s │ 0.002GB │        42 │      1.0 │
  └───────────┴──────────┴─────────┴───────────┴──────────┘

  Errors:

  Trotter   vs RK4 rel=2.85e-3, max=2.32e-1
  TrotterTS vs RK4 rel=3.96e-3, max=2.11e-1

  Caution: N=64 TrotterTS is extremely fast because it retained only 42 terms,
  so this is a very aggressive component truncation regime.

  ---
  Summary

  - TrotterTS is much faster for larger N when componenttol is modest (0.999)
  because it keeps far fewer terms.
  - Accuracy is generally comparable to expanded Trotter at coarse tolerances,
  but often slightly worse at small exact cases.
  - At tighter settings (dt=0.05, tol=0.9999), TrotterTS can be close to
  expanded Trotter accuracy and faster for N=32, but not for small N.
  - TrotterTS allocations are low when component truncation is aggressive, but
  can exceed RK4 for tighter/small cases.

  The clearest promising regime from this sweep is:

  N=32, dt=0.05, tol=0.9999:
  TrotterTS ~5.3x faster than RK4
  rel vs RK4 ~1.5e-3
  final len 503 vs 4096

✻ Baked for 2m 41s

※ recap: We’re developing a TS-native Trotter method and currently comparing TrotterTS against RK4 and expanded Trotter across parameters. Next, decide whether the current component-truncation accuracy/speed tradeoff is acceptable or adjust the method. (disable recaps in /config)

❯ I see. So it seems the drastic speedup we get by freezing the tail components 
  leads to poor accuracy. Is there a method that applies the update to the      
  tail orbits in its entirity? That is, not individually to each orbit, but     
  rather to the entire tail in one go                                           

● Yes. That is a very natural next idea.

  Right now, for one orbit flow we split active components into:

  head components: exact component exponential
  tail components: frozen

  The problem is that freezing the tail is dynamically wrong.

  Instead, we could treat the tail as one aggregate subspace and apply a cheap
  approximate update to the entire tail at once.

  That would avoid:

  one expensive exact exponential per tail component

  while still evolving the tail instead of freezing it.

  Conceptual decomposition

  Let the active operator be decomposed into components under one orbit
  Liouvillian:

  O = O_H + O_T,

  where:

  - O_H = sum of high-weight components,
  - O_T = sum of tail components.

  Because components are invariant under the orbit Liouvillian:

  L_a(C_i) \subseteq C_i,

  we also have:

  L_a(O_T)

  stays in the union of tail components.

  So we can update the tail as one operator:

  O_T' \approx \Phi_{\mathrm{cheap}}(\tau L_a) O_T.

  Then:

  O' = \sum_{C \in H} e^{\tau L_{a,C}} O_C + O_T'.

  Cheap tail update options

  1. One global Taylor step on the tail

  O_T' = O_T + \tau L_a O_T.

  or second order:

  O_T' = O_T + \tau L_a O_T + \frac{\tau^2}{2} L_a^2 O_T.

  This uses one or two calls to:

  orbit_liouvillian(Ha, O_tail)

  rather than many component matrix exponentials.

  This is likely the simplest.

  It is still not exact Trotter for the tail, but it is much better than
  freezing.

  2. One small Krylov update on the entire tail

  Instead of Krylov per component, do:

  O_T' \approx e^{\tau L_a} O_T

  with Krylov dimension 1–2 or 3.

  This applies the orbit flow approximately to the entire tail in one go.

  Cost:

  k orbit_liouvillian calls on O_tail

  rather than:

  #tail components dense exponentials

  3. RK-like full tail update

  For the tail only:

  use RK2 / midpoint / Heun under L_a

  But this is basically a polynomial tail subflow.

  4. Freeze plus correction

  A hybrid:

  O_T' = O_T + \tau L_a O_T

  but then trim aggressively.

  This may be enough.

  Why this may work

  Tail components collectively have small coefficient norm, so they do not need
  exact high-order evolution.

  But freezing them creates a systematic phase/dynamics error. A first- or
  second-order aggregate update may remove most of that error while preserving
  the component speedup.

  The key is:

  update the tail once as an OperatorTS,
  not component-by-component.

  Implementation shape

  Inside _orbit_flow:

  1. Build component data and sort by weight.
  2. For components in the head:
    - exact component exponential.
  3. For components in the tail:
    - accumulate their active coefficients into one Otail.
  4. Apply cheap tail update:

  Ot = Otail + dt * orbit_liouvillian(Ha, Otail)

  4. or Taylor-2:

  K1 = orbit_liouvillian(Ha, Otail)
  K2 = orbit_liouvillian(Ha, K1)
  Ot = Otail + dt*K1 + dt^2/2*K2
  5. Return:

  Ohead + Ot

  Why this is different from the previous tail1 attempt

  Previously, I updated each low-weight component separately:

  for every tail component:
      coeffs += dt * A_component * coeffs

  That still touched all components and was expensive.

  This new version would do:

  Otail = sum of all tail components
  Otail += dt * L_a(Otail)

  one global tail update.

  So the number of tail components no longer directly controls tail update cost.

  What to test

  Modes:

  tail = :freeze
  tail = :euler
  tail = :taylor2
  tail = :krylov2

  At fixed:

  componenttol = 0.999 or 0.9999

  Measure:

  - accuracy vs exact/RK4,
  - time,
  - allocations,
  - final length.

  Expected:

  - :freeze fastest, less accurate.
  - :euler slightly slower, hopefully more accurate.
  - :taylor2 slower but maybe close to full exact component update.
  - :krylov2 maybe best accuracy/cost compromise.

  This is worth testing.

✻ Baked for 34s

❯ before we do this, is it worth it to apply component batching for the orbits, 
  and then checking for weights of the batches?                                 

● Yes, it is worth considering — and it might be better than weighting
  individual components.

  The idea would be:

  group components by shared component type / signature
  then decide which batches to update exactly based on total batch weight

  instead of:

  rank every individual component by its own weight

  This could reduce overhead and maybe improve accuracy.

  Why batching may help

  Right now we may have hundreds of components, many of which have the same
  structure.

  Example from the signature check:

  N=32, orbit 3:
    697 components
    44 spectral signatures
    top signature appears 398 times

  That suggests many components are structurally identical.

  If we batch them, we can:

  1. Build one component matrix per signature.
  2. Compute one exponential per signature and timestep.
  3. Apply it to many coefficient vectors.

  This reduces the overhead of:

  many small exp(A)

  and may make exact updates cheaper.

  Batch-level weighting

  Instead of component weights:

  w_C = \|c_C\|^2

  define batch weight:

  w_B = \sum_{C \in B} w_C.

  Then keep/update batches by descending total weight.

  This can be better because many individually small components of the same type
  may collectively matter.

  For example:

  one component weight: small
  500 similar components: total weight significant

  Component-wise truncation might freeze all of them individually. Batch-wise
  truncation would keep the whole family.

  That could improve accuracy.

  But there is a risk

  Batching by structural type is useful only if applying the same matrix
  exponential works across all components.

  For that, we need a consistent local ordering. If two components have the same
  graph structure but different node labels, then:

  A_component_1
  A_component_2

  must be identical up to permutation.

  A robust implementation needs a canonical component labeling/signature, not
  just same eigenvalues.

  Spectral signatures can collide: two non-isomorphic matrices can have the same
  eigenvalues.

  Degree signatures are even weaker.

  So for correctness, batching requires a canonical graph/matrix form.

  Minimal experiment

  Before implementing batching, measure batch weights by approximate signatures.

  Use two signatures:

  1. Exact matrix signature in current component order, rounded:
    - this only groups components whose matrices happen to have identical
  ordering.
    - cheap and safe for sharing exp(A) only if basis ordering matches.
  2. Spectral signature, for diagnostic only:
    - tells us potential batching upper bound.
    - not safe alone for applying updates.

  For each orbit and state, compute:

  number of components
  number of exact-matrix signatures
  number of spectral signatures
  weight concentration by component
  weight concentration by signature batch

  Then compare:

  components needed for 99.9% weight
  batches needed for 99.9% weight

  If batch weights are more concentrated, batching is promising.

  Potential outcomes

  Good outcome

  700 components
  50 exact matrix signatures
  top 5 signatures carry 99% weight

  Then batching is very promising.

  Bad outcome

  700 components
  600 exact matrix signatures
  batch weights similar to component weights

  Then batching is not worth it.

  Middle outcome

  few spectral signatures but many exact matrix signatures

  Then batching might require canonical relabeling, which is more involved.

  Should we do it before tail update?

  Yes. It is a cleaner structural optimization:

  - It preserves exact component updates for retained batches.
  - It may reduce the number of exponentials and bookkeeping operations.
  - It may reduce the need for freezing/tail approximation.

  So I would run the batch-weight diagnostic first. That tells us whether
  component batching is a viable path.

✻ Worked for 32s

❯ go                                                                            

● Write(~/.claude/jobs/513a65b3/tmp/batch_weight_diagnostic.jl)
  ⎿  Wrote 128 lines to
     ../../../.claude/jobs/513a65b3/tmp/batch_weight_diagnostic.jl
       1 using PauliStrings
       2 using LinearAlgebra
       3 using Printf
       5 
       6 function MFIM(N,h)
       7     H=Operator(N); H += -h,"X",1; H += -h/2,"Z",2; H += "Z",1,"Z
         ",2; OperatorTS{(N,)}(H)
       8 end
         ; end
      10 function grow(H,O0,M,steps)
     … +149 lines (ctrl+o to expand)

  Ran 1 shell command (ctrl+o to expand)

● Ran the batch-weight diagnostic.

  Script:

  /home/dhruva/.claude/jobs/513a65b3/tmp/batch_weight_diagnostic.jl

  It compares:

  - individual component weight selection,
  - exact-matrix-signature batch selection,
  - spectral-signature batch selection,
  - degree-signature batch selection.

  Key finding

  Batching by component signature could reduce the number of exact component 
  exponentials substantially.

  For the same retained weight, we often need far fewer batches than individual
  components.

  ---
  N=16, M=10

  Orbit 1

  217 components
  25 exact matrix signatures

  Components/batches needed:

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │               5 / 217 │               2 / 25 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │              12 / 217 │               2 / 25 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │              19 / 217 │               2 / 25 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │              30 / 217 │               3 / 25 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 2

  229 components
  48 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     5 │                    3 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                    13 │                    5 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    24 │                    6 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    38 │                    9 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 3

  251 components
  42 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     4 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                     8 │                    4 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    19 │                    6 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    36 │                    7 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  ---
  N=32, M=12

  Orbit 1

  670 components
  134 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     5 │                    3 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                    13 │                    5 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    24 │                    6 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    38 │                    9 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 2

  782 components
  123 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     4 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                     8 │                    3 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    19 │                    6 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    36 │                    7 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 3

  697 components
  45 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     5 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                    12 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    19 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    30 │                    3 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  This is especially promising: orbit 3 needs only 2–3 exact matrix signatures
  to cover almost all weight.

  ---
  N=64, M=14

  Orbit 1

  1869 components
  496 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     5 │                    3 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                    13 │                    4 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    24 │                    5 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    38 │                    8 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 2

  2466 components
  545 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     4 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                     8 │                    3 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    19 │                    5 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    36 │                    7 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  Orbit 3

  2168 components
  157 exact matrix signatures

  ┌─────────────────┬───────────────────────┬──────────────────────┐
  │ retained weight │ individual components │ exact matrix batches │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             90% │                     5 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │             99% │                    12 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │           99.9% │                    19 │                    2 │
  ├─────────────────┼───────────────────────┼──────────────────────┤
  │          99.99% │                    30 │                    3 │
  └─────────────────┴───────────────────────┴──────────────────────┘

  ---
  Interpretation

  Batching looks worthwhile.

  Component-wise truncation says:

  keep top ~30–40 individual components for 99.99%

  Batch-wise truncation says:

  keep top ~3–9 exact matrix signatures for 99.99%

  This matters because if all components in one batch share the same component
  matrix, then:

  one exp(A)
  many coefficient vectors

  instead of:

  many exp(A_i)

  Also, batch weighting may be more accurate than component weighting because
  many individually small but structurally identical components can collectively
  matter.

  ---
  Caveat

  The “exact matrix signature” used here groups components whose matrix entries
  match in their current component ordering. This is safe for sharing the matrix
  exponential only if the ordering is consistent. It is conservative: it may
  miss isomorphic components with different ordering.

  So this is a lower bound on batching potential. Spectral/degree signatures
  show even more possible grouping, but they are not safe for direct reuse
  without canonical relabeling.

  ---
  Next implementation idea

  Instead of selecting individual components by weight:

  sort components by component weight

  select batches by exact matrix signature:

  sort signatures by total batch weight

  Then for every component whose signature is retained:

  apply exact component exponential

  For components in unretained signatures:

  freeze

  This should reduce the number of exponentials while retaining whole component
  families, potentially improving both speed and accuracy.

✻ Baked for 1m 32s

❯ I see, but this all depends on the ability to apply the exponent to the       
  entire batch at once?                                                         

● Yes — batching is only truly useful if we can reuse the same exponential
  action across all components in the batch.

  But there are two levels of “apply to the whole batch at once.”

  Level 1: reuse the same exp(A)

  Suppose many components have the same local Liouvillian matrix:

  A_C = A_{\mathrm{sig}}

  Then we compute:

  E = exp(dt * A_sig)

  once, and for each component in the batch do:

  coeffs_out = E * coeffs_in

  This already saves the expensive dense exponential construction.

  It does not eliminate per-component gather/scatter or matvec costs.

  So cost becomes:

  one exp per signature
  + one small matvec per component
  + gather/scatter per component

  This is useful if exp(A) dominates.

  Level 2: batched matvec for all components

  If many components share the same matrix E, stack their coefficient vectors as
  columns:

  C = [c_1 \; c_2 \; c_3 \; \cdots]

  Then apply:

  C' = E C.

  In Julia/BLAS terms:

  Cout = E * Cin

  where Cin is a dense matrix with one column per component in the batch.

  This applies the exponential action to the whole batch in one dense
  matrix-matrix multiply.

  That is much better than many tiny matrix-vector multiplies.

  So batching is most useful if we do Level 2.

  What is required for this to be valid

  Components in a batch must share:

  1. same component size,
  2. same local Liouvillian matrix in the same canonical ordering,
  3. same dt.

  Then a single E applies to all coefficient columns.

  Our “exact matrix signature” checks for same matrix in current ordering, so it
  is safe for this kind of batching.

  If we later want to batch isomorphic components with permuted orderings, we
  need canonical relabeling. That is harder.

  Cost model

  Without batching:

  for each component:
      build/lookup E_i
      coeff_out_i = E_i * coeff_in_i

  With Level 1 batching:

  for each signature:
      build/lookup E_sig once
      for each component in batch:
          coeff_out_i = E_sig * coeff_in_i

  With Level 2 batching:

  for each signature:
      build/lookup E_sig once
      Cin = matrix of all coeff vectors in batch
      Cout = E_sig * Cin

  The last one is where real speedups may appear.

  Does batching help if we do not use Level 2?

  Maybe, but less.

  Our current component code already caches each component’s exp(A) plan. If
  many components have the same matrix but are stored as separate plans, then we
  are still computing duplicate exponentials. Level 1 avoids that.

  But if component sizes are very small, the bigger overhead may be
  gather/scatter, not exp(A). Then Level 2 is needed.

  Why batch weights matter

  If we select by batch weight, we update whole families:

  all components sharing signature S

  instead of just a few individual components.

  That could improve accuracy because it avoids freezing many low-weight but
  collectively important components.

  But it could also increase work if a large batch has many components. That is
  why batched matvec matters.

  Short answer

  Yes. Batching is only a serious optimization if we can apply:

  E_{\mathrm{sig}}

  to all components in a signature together.

  The ideal implementation is:

  Cout = E_sig * Cin

  where each column of Cin is one component’s coefficient vector.

✻ Cogitated for 32s

※ recap: We’re trying to make TS Trotter faster by using translation-orbit component updates instead of expanded gates. Next, test batching identical component matrices so one exponential can update many components at once. (disable recaps in /config)

❯ are the components we have in a batch already of the same size with the       
  component similarity method you already chose?                                

● Yes, with the exact matrix signature I used in the diagnostic, components in a
  batch are already the same size.

  The signature was effectively:

  (size(A, 1), nonzero matrix entries in current ordering)

  So two components are in the same exact-matrix batch only if:

  1. they have the same size(A,1),
  2. they have the same nonzero matrix pattern,
  3. they have the same matrix values,
  4. in the same local ordering.

  Therefore, they can share the same dense exponential:

  E = exp(dt * A)

  and they can be batched as columns:

  Cin  :: Matrix{ComplexF64}  # size n_component × n_batch
  Cout = E * Cin

  where every column corresponds to one component in that batch.

  Caveat

  This is true for the exact matrix signature, not for the weaker signatures.

  - Exact matrix signature: same size and same matrix in same ordering — safe
  for batched application.
  - Spectral signature: same eigenvalues only — not safe.
  - Degree signature: same degree pattern only — not safe.

  So if we use the exact matrix signature, yes, batching is directly valid.

  The limitation

  The exact matrix signature is conservative. It may fail to group components
  that are actually isomorphic but ordered differently.

  But that is fine for a first implementation: safe batching first, canonical
  relabeling later if needed.
✻ Cogitated for 22s

❯ Ok then go ahead and implement using exact matrix signatures.                 

● API Error: Request rejected (429) · Sie haben $20.496685 / $20.0 
  verbraucht und damit ihr monatliches Budget erreicht.
  Normalerweise wird das Budget zu Beginn jedes Monats zurückgesetzt.
  Mehr Informationen über das Thema Budgetierung finden Sie im 
  Hilfe-Wiki unter 
  https://hilfe.uni-paderborn.de/Budgets_im_AI-Chat_und_AI-Gateway

✻ Sautéed for 3m 15s

● Auto mode lets Claude handle permission prompts automatically — Claude
  checks each tool call for risky actions and prompt injection before 
  executing. Actions Claude identifies as safe are executed, while 
  actions Claude identifies as risky are blocked and Claude may try a 
  different approach. Ideal for long-running tasks. Sessions are 
  slightly more expensive. Claude can make mistakes that allow harmful 
  commands to run, it's recommended to only use in isolated 
  environments. Shift+Tab to change mode.