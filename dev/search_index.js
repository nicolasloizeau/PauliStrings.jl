var documenterSearchIndex = {"docs":
[{"location":"","page":"Getting started","title":"Getting started","text":"(Image: Build Status) (Image: )","category":"page"},{"location":"#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"PauliStrings.jl is a Julia package for many-body quantum mechanics with Pauli string represented as binary integers (as in https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318).","category":"page"},{"location":"#Installation","page":"Getting started","title":"Installation","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"using Pkg; Pkg.add(url=\"https://github.com/nicolasloizeau/PauliStrings.jl\") or ] add https://github.com/nicolasloizeau/PauliStrings.jl","category":"page"},{"location":"#Initializing-an-operator","page":"Getting started","title":"Initializing an operator","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"Import the library and initialize a operator of 4 qubits","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"using PauliStrings\nimport PauliStrings as ps\nH = ps.Operator(4)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Add a Pauli strings to the operator","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H += \"XYZ1\"\nH += \"1YZY\"","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"julia> H\n(1.0 - 0.0im) XYZ1\n(1.0 - 0.0im) 1YZY","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Add a Pauli string with a coeficient","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H += -1.2,\"XXXZ\" #coeficient can be complex","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Add a 2-qubit string coupling qubits i and j with X and Y:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H += 2, \"X\", i, \"Y\", j # with a coeficient=2\nH += \"X\", i, \"Y\", j # with a coeficient=1","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Add a 1-qubit string:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H += 2, \"Z\", i # with a coeficient=2\nH += \"Z\", i # with a coeficient=1\nH += \"S+\", i","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Supported sites operators are X, Y, Z, Sx=X2, Sy=Y2, Sz=Z2, S+=(X+iY)2, S-=(X-iY)2.","category":"page"},{"location":"#Basic-Algebra","page":"Getting started","title":"Basic Algebra","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"The Operator type supports the +,-,* operators with other Operators and Numbers:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H3 = H1*H2\nH3 = H1+H2\nH3 = H1-H2\nH3 = H1+2 # adding a scalar is equivalent to adding the unit times the scalar\nH = 5*H # multiply operator by a scalar","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Trace : ps.trace(H)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Frobenius norm : ps.opnorm(H)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Conjugate transpose : ps.dagger(H)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Number of terms: length(H)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Commutator: ps.com(H1, H2). This is much faster than H1*H2-H2*H1","category":"page"},{"location":"#Print-and-export","page":"Getting started","title":"Print and export","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"print shows a list of terms with coeficients e.g :","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"julia> println(H)\n(10.0 - 0.0im) 1ZZ\n(5.0 - 0.0im) 1Z1\n(15.0 + 0.0im) XYZ\n(5.0 + 0.0im) 1YY","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Export a list of strings with coeficients:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"coefs, strings = ps.op_to_strings(H)","category":"page"},{"location":"#Truncate,-Cutoff,-Trim,-Noise","page":"Getting started","title":"Truncate, Cutoff, Trim, Noise","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"ps.truncate(H,M) removes Pauli strings longer than M (returns a new Operator) ps.cutoff(H,c) removes Pauli strings with coeficient smaller than c in absolute value (returns a new Operator) ps.trim(H,N) keeps the first N trings with higest weight (returns a new Operator) ps.prune(H,alpha) keeps terms with probability 1-exp(-alpha*abs(c)) (returns a new Operator)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"ps.add_noise(H,g) adds depolarizing noise that make each strings decay like e^gw where w is the lenght of the string. This is usefull when used with trim to keep the number of strings manageable during time evolution.","category":"page"},{"location":"#Time-evolution","page":"Getting started","title":"Time evolution","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"ps.rk4(H, O, dt; hbar=1, heisenberg=false) performs a step of Runge Kutta and returns the new updated O(t+dt)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H can be an Operator, or a function that takes a time and return an Operator. In case H is a function, a time also needs to be passed to rk4(H, O, dt, t). O is an Observable or a density matrix to time evolve. If evolving an observable in the heisenberg picture, set heisenberg=true.","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"An example is in time_evolve_example.jl. The following will time evolve O in the Heisenberg picture. At each step, we add depolarizing noise and trim the operator to keep the number of strings manageable","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"function evolve(H, O, M, times, noise)\n    dt = times[2]-times[1]\n    for t in times\n        O = ps.rk4(H, O, dt; heisenberg=true, M=M) #preform one step of rk4, keep only M strings\n        O = ps.add_noise(O, noise*dt) #add depolarizingn noise\n        O = ps.trim(O, M) # keep the M strings with the largest weight\n    end\n    return O\nend","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Time evolution of the spin correlation function textupTr(Z_1(0)Z_1(t)) in the chaotic spin chain. Check timeevolveexample.jl to reproduce the plot. (Image: plot)","category":"page"},{"location":"#Lanczos","page":"Getting started","title":"Lanczos","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"Compute lanczos coeficients","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"bs = ps.lanczos(H, O, steps, nterms)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"H : Hamiltonian","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"O : starting operator","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"nterms : maximum number of terms in the operator. Used by trim at every step","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Results for X in XX from https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 :","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"(Image: plot)","category":"page"},{"location":"documentation/#Documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"documentation/#Basics","page":"Documentation","title":"Basics","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Operator(N::Int)","category":"page"},{"location":"documentation/#PauliStrings.Operator-Tuple{Int64}","page":"Documentation","title":"PauliStrings.Operator","text":"Operator(N::Int)\n\nInitialize an empty operator on N qubits\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Base.length(o::Operator)","category":"page"},{"location":"documentation/#Base.length-Tuple{Operator}","page":"Documentation","title":"Base.length","text":"Base.length(o::Operator)\n\nNumber of pauli strings in an operator\n\nExample\n\njulia> A = Operator(4)\njulia> A += \"X111\"\njulia> A += \"XYZ1\"\njulia> A += 2, \"Y\", 4\njulia> length(A)\n3\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"eye(N::Int)","category":"page"},{"location":"documentation/#PauliStrings.eye-Tuple{Int64}","page":"Documentation","title":"PauliStrings.eye","text":"eye(N::Int)\n\nIdentity operator on N qubits\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Base.:+(o::Operator, args::Tuple{Number, Vararg{Any}})\nBase.:+(o::Operator, args::Tuple{Vararg{Any}})\nBase.:+(o::Operator, term::Tuple{Number, String})\nBase.:+(o::Operator, term::String)","category":"page"},{"location":"documentation/#Base.:+-Tuple{Operator, Tuple{Number, Vararg{Any, N} where N}}","page":"Documentation","title":"Base.:+","text":"Base.:+(o::Operator, args::Tuple{Number, Vararg{Any}})\nBase.:+(o::Operator, args::Tuple{Vararg{Any}})\nBase.:+(o::Operator, term::Tuple{Number, String})\nBase.:+(o::Operator, term::String)\n\nMain functions to contruct spin operators. Identical signatures are available for -.\n\nExamples\n\nk-local terms can be added by adding a tuple to the operator. The first element of the tuple is an optional coeficient. The other element are couples (symbol,site) where symbol can be \"X\", \"Y\", \"Z\", \"Sx\", \"Sy\", \"Sz\", \"S+\", \"S-\" and site is an integer specifying the site on wich the symbol is acting.\n\nA = Operator(4)\nA += 2, \"X\",1,\"X\",2\nA += 3, \"Y\",1,\"X\",2\nA += \"X\",3,\"X\",4\nA += 4,\"Z\",3\nA += 5.2,\"X\",1,\"Y\",2,\"Z\",3\n\njulia> A\n(4.0 + 0.0im) 11Z1\n(3.0 - 0.0im) YX11\n(1.0 + 0.0im) 11XX\n(2.0 + 0.0im) XX11\n(5.2 - 0.0im) XYZ1\n\nFull strings can also be added:\n\nA = Operator(4)\nA += 2, \"1XXY\"\nA += 2im, \"11Z1\"\n\njulia> A\n(0.0 + 2.0im) 11Z1\n(2.0 - 0.0im) 1XXY\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Operator, Tuple}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Operator, Tuple{Number, String}}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\nBase.:+(o::Operator, args::Tuple{Number, Vararg{Any}})\nBase.:+(o::Operator, args::Tuple{Vararg{Any}})\nBase.:+(o::Operator, term::Tuple{Number, String})\nBase.:+(o::Operator, term::String)\n\nMain functions to contruct spin operators. Identical signatures are available for -.\n\nExamples\n\nk-local terms can be added by adding a tuple to the operator. The first element of the tuple is an optional coeficient. The other element are couples (symbol,site) where symbol can be \"X\", \"Y\", \"Z\", \"Sx\", \"Sy\", \"Sz\", \"S+\", \"S-\" and site is an integer specifying the site on wich the symbol is acting.\n\nA = Operator(4)\nA += 2, \"X\",1,\"X\",2\nA += 3, \"Y\",1,\"X\",2\nA += \"X\",3,\"X\",4\nA += 4,\"Z\",3\nA += 5.2,\"X\",1,\"Y\",2,\"Z\",3\n\njulia> A\n(4.0 + 0.0im) 11Z1\n(3.0 - 0.0im) YX11\n(1.0 + 0.0im) 11XX\n(2.0 + 0.0im) XX11\n(5.2 - 0.0im) XYZ1\n\nFull strings can also be added:\n\nA = Operator(4)\nA += 2, \"1XXY\"\nA += 2im, \"11Z1\"\n\njulia> A\n(0.0 + 2.0im) 11Z1\n(2.0 - 0.0im) 1XXY\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Operator, String}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Operations","page":"Documentation","title":"Operations","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"add(o1::Operator, o2::Operator)\nBase.:+(o1::Operator, o2::Operator)\nBase.:+(o::Operator, a::Number)\nBase.:+(a::Number, o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.add-Tuple{Operator, Operator}","page":"Documentation","title":"PauliStrings.add","text":"add(o1::Operator, o2::Operator)\nBase.:+(o1::Operator, o2::Operator)\nBase.:+(o::Operator, a::Number)\nBase.:+(a::Number, o::Operator)\n\nAdd two operators together or add a number to an operator\n\nExample\n\nA = Operator(4)\nA += \"XYZ1\"\nA += 1, \"Y\", 4\nB = Operator(4)\nB += 2, \"Y\", 2, \"Y\", 4\nB += 1, \"Z\", 3\n\njulia> A\n(1.0 - 0.0im) 111Y\n(1.0 - 0.0im) XYZ1\n\njulia> B\n(1.0 + 0.0im) 11Z1\n(2.0 - 0.0im) 1Y1Y\n\njulia> A+B\n(1.0 + 0.0im) 11Z1\n(2.0 - 0.0im) 1Y1Y\n(1.0 - 0.0im) 111Y\n(1.0 - 0.0im) XYZ1\n\njulia> A+5\n(1.0 - 0.0im) 111Y\n(1.0 - 0.0im) XYZ1\n(5.0 + 0.0im) 1111\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Operator, Operator}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Operator, Number}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:+-Tuple{Number, Operator}","page":"Documentation","title":"Base.:+","text":"+(x, y...)\n\nAddition operator. x+y+z+... calls this function with all arguments, i.e. +(x, y, z, ...).\n\nExamples\n\njulia> 1 + 20 + 4\n25\n\njulia> +(1, 20, 4)\n25\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Base.:*(o1::Operator, o2::Operator)\nBase.:*(o::Operator, a::Number)\nBase.:*(a::Number, o::Operator)","category":"page"},{"location":"documentation/#Base.:*-Tuple{Operator, Operator}","page":"Documentation","title":"Base.:*","text":"Base.:*(o1::Operator, o2::Operator)\nBase.:*(o::Operator, a::Number)\nBase.:*(a::Number, o::Operator)\n\nMultiply two operators together or an operator with a number\n\nExample\n\nA = Operator(4)\nA += \"XYZ1\"\nA += 1, \"Y\", 4\nB = Operator(4)\nB += 2, \"Y\", 2, \"Y\", 4\nB += 1, \"Z\", 3\n\njulia> A\n(1.0 - 0.0im) 111Y\n(1.0 - 0.0im) XYZ1\n\n\njulia> B\n(1.0 + 0.0im) 11Z1\n(2.0 - 0.0im) 1Y1Y\n\njulia> A*B\n(2.0 - 0.0im) X1ZY\n(1.0 - 0.0im) 11ZY\n(2.0 - 0.0im) 1Y11\n(1.0 - 0.0im) XY11\n\njulia> A*5\n(5.0 - 0.0im) 111Y\n(5.0 - 0.0im) XYZ1\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:*-Tuple{Operator, Number}","page":"Documentation","title":"Base.:*","text":"*(x, y...)\n\nMultiplication operator. x*y*z*... calls this function with all arguments, i.e. *(x, y, z, ...).\n\nExamples\n\njulia> 2 * 7 * 8\n112\n\njulia> *(2, 7, 8)\n112\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:*-Tuple{Number, Operator}","page":"Documentation","title":"Base.:*","text":"*(x, y...)\n\nMultiplication operator. x*y*z*... calls this function with all arguments, i.e. *(x, y, z, ...).\n\nExamples\n\njulia> 2 * 7 * 8\n112\n\njulia> *(2, 7, 8)\n112\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Base.:-(o::Operator)\nBase.:-(o1::Operator, o2::Operator)\nBase.:-(o::Operator, a::Real)\nBase.:-(a::Real, o::Operator)","category":"page"},{"location":"documentation/#Base.:--Tuple{Operator}","page":"Documentation","title":"Base.:-","text":"Base.:-(o::Operator)\nBase.:-(o1::Operator, o2::Operator)\nBase.:-(o::Operator, a::Real)\nBase.:-(a::Real, o::Operator)\n\nSubtraction between operators and numbers\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:--Tuple{Operator, Operator}","page":"Documentation","title":"Base.:-","text":"-(x, y)\n\nSubtraction operator.\n\nExamples\n\njulia> 2 - 3\n-1\n\njulia> -(2, 4.5)\n-2.5\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:--Tuple{Operator, Real}","page":"Documentation","title":"Base.:-","text":"-(x, y)\n\nSubtraction operator.\n\nExamples\n\njulia> 2 - 3\n-1\n\njulia> -(2, 4.5)\n-2.5\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Base.:--Tuple{Real, Operator}","page":"Documentation","title":"Base.:-","text":"-(x, y)\n\nSubtraction operator.\n\nExamples\n\njulia> 2 - 3\n-1\n\njulia> -(2, 4.5)\n-2.5\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000)","category":"page"},{"location":"documentation/#PauliStrings.com-Tuple{Operator, Operator}","page":"Documentation","title":"PauliStrings.com","text":"com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000)\n\nCommutator of two operators\n\nExample\n\njulia> A = Operator(4)\njulia> A += \"X111\"\njulia> B = Operator(4)\njulia> B += \"Z111\"\njulia> B += \"XYZ1\"\njulia> com(A,B)\n(0.0 - 2.0im) Y111\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"trace(o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.trace-Tuple{Operator}","page":"Documentation","title":"PauliStrings.trace","text":"trace(o::Operator)\n\nTrace of an operator\n\nExample\n\njulia> A = Operator(4)\njulia> A += 2,\"1111\"\njulia> A += 3,\"XYZ1\"\njulia> trace(A)\n32.0 + 0.0im\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"opnorm(o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.opnorm-Tuple{Operator}","page":"Documentation","title":"PauliStrings.opnorm","text":"opnorm(o::Operator)\n\nFrobenius norm\n\nExample\n\njulia> A = Operator(4)\njulia> A += 2,\"X\",2\njulia> A += 1,\"Z\",1,\"Z\",3\njulia> opnorm(A)\n8.94427190999916\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"dagger(o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.dagger-Tuple{Operator}","page":"Documentation","title":"PauliStrings.dagger","text":"dagger(o::Operator)\n\nConjugate transpose\n\nExample\n\nA = Operator(3)\nA += 1im,\"X\",2\nA += 1,\"Z\",1,\"Z\",3\n\njulia> A\n\n(1.0 + 0.0im) Z1Z\n(0.0 + 1.0im) 1X1\n\n\njulia> dagger(A)\n(1.0 - 0.0im) Z1Z\n(0.0 - 1.0im) 1X1\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"ptrace(o::Operator, keep::Vector{Int})","category":"page"},{"location":"documentation/#PauliStrings.ptrace-Tuple{Operator, Vector{Int64}}","page":"Documentation","title":"PauliStrings.ptrace","text":"ptrace(o::Operator, keep::Vector{Int})\n\nPartial trace.\n\nkeep is list of qubits indices to keep starting at 1 note that this still returns an operator of size N and doesnt permute the qubits this only gets rid of Pauli strings that have no support on keep and add their coeficient*2^N to the identity string\n\nExample\n\nA = Operator(5)\nA += \"XY1XZ\"\nA += \"XY11Z\"\n\njulia> ptrace(A, [3,4])\n(1.0 - 0.0im) XY1XZ\n(8.0 - 0.0im) 11111\n\njulia> ptrace(A, [1,5])\n(1.0 - 0.0im) XY1XZ\n(1.0 - 0.0im) XY11Z\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Random-operators","page":"Documentation","title":"Random operators","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"rand_local1(N::Int)","category":"page"},{"location":"documentation/#PauliStrings.rand_local1-Tuple{Int64}","page":"Documentation","title":"PauliStrings.rand_local1","text":"rand_local2(N::Int)\n\nRandom 1-local operator\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"rand_local2(N::Int)","category":"page"},{"location":"documentation/#PauliStrings.rand_local2-Tuple{Int64}","page":"Documentation","title":"PauliStrings.rand_local2","text":"rand_local2(N::Int)\n\nRandom 2-local operator\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Truncation-and-noise","page":"Documentation","title":"Truncation and noise","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"truncate(o::Operator, N::Int; keepnorm::Bool = false)","category":"page"},{"location":"documentation/#Base.truncate-Tuple{Operator, Int64}","page":"Documentation","title":"Base.truncate","text":"truncate(o::Operator, N::Int; keepnorm::Bool = false)\n\nRemove all terms of length > N. Keep all terms of length <= N. i.e remove all M-local terms with M>N\n\nExample\n\nA = Operator(4)\nA += \"X\",1,\"X\",2\nA += \"Z\",1,\"Z\",2,\"Z\",4\n\njulia> A\n(1.0 + 0.0im) ZZ1Z\n(1.0 + 0.0im) XX11\n\njulia> ps.truncate(A,2)\n(1.0 + 0.0im) XX11\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"trim(o::Operator, N::Int; keepnorm::Bool = false, keep::Operator=Operator(N))","category":"page"},{"location":"documentation/#PauliStrings.trim-Tuple{Operator, Int64}","page":"Documentation","title":"PauliStrings.trim","text":"trim(o::Operator, N::Int; keepnorm::Bool = false, keep::Operator=Operator(N))\n\nKeep the first N terms with largest coeficients.\n\nkeepnorm is set to true to keep the norm of o.\n\nkeep is an operator that specify a set of strings that cannot be removed\n\nExample\n\nA = Operator(4)\nA += 1,\"XXXX\"\nA += 2,\"XX11\"\nA += 3,\"XX1X\"\nA += 4,\"ZZXX\"\nB = Operator(4)\nB += 1,\"XX11\"\nB += 1,\"XX1X\"\n\njulia> trim(A,2)\n(4.0 + 0.0im) ZZXX\n(3.0 + 0.0im) XX1X\n\njulia> trim(A,2;keep=B)\n(4.0 + 0.0im) ZZXX\n(3.0 + 0.0im) XX1X\n(2.0 + 0.0im) XX11\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"prune(o::Operator, alpha::Real; keepnorm::Bool = false)","category":"page"},{"location":"documentation/#PauliStrings.prune-Tuple{Operator, Real}","page":"Documentation","title":"PauliStrings.prune","text":" prune(o::Operator, alpha::Real; keepnorm::Bool = false)\n\nKeep terms with probability 1-exp(-alpha*abs(c)) where c is the weight of the term\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"cutoff(o::Operator, epsilon::Real; keepnorm::Bool = false)","category":"page"},{"location":"documentation/#PauliStrings.cutoff-Tuple{Operator, Real}","page":"Documentation","title":"PauliStrings.cutoff","text":"cutoff(o::Operator, epsilon::Real; keepnorm::Bool = false)\n\nRemove all terms with weight < epsilon\n\nExample\n\nA = Operator(4)\nA += 1,\"XXXX\"\nA += 2,\"XX11\"\nA += 3,\"XX1X\"\nA += 4,\"ZZXX\"\n\njulia> cutoff(A, 2.5)\n(3.0 + 0.0im) XX1X\n(4.0 + 0.0im) ZZXX\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"add_noise(o::Operator, g::Real)","category":"page"},{"location":"documentation/#PauliStrings.add_noise-Tuple{Operator, Real}","page":"Documentation","title":"PauliStrings.add_noise","text":"add_noise(o::Operator, g::Real)\n\nAdd depolarizing noise that make the long string decays. g is the noise amplitude.\n\nExample\n\nA = add_noise(A, 0.1)\n\nReference\n\nhttps://arxiv.org/pdf/2407.12768\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Algorithms","page":"Documentation","title":"Algorithms","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)","category":"page"},{"location":"documentation/#PauliStrings.lanczos-Tuple{Operator, Operator, Int64, Int64}","page":"Documentation","title":"PauliStrings.lanczos","text":"lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)\n\nComputer the first steps lanczos coeficients for Hamiltonian H and initial operator O\n\nAt every step, the operator is trimed with PauliStrings.trim and only nterms are kept.\n\nIf H and O are 1D-translation-invariant, it is possible to provide a single local term of O and set localop=true\n\nUsing maxlength speeds up the commutator by only keeping terms of length <= maxlength\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))","category":"page"},{"location":"documentation/#PauliStrings.rk4-Tuple{Operator, Operator, Real}","page":"Documentation","title":"PauliStrings.rk4","text":"rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))\n\nSingle step of Runge–Kutta-4 with time independant Hamiltonian. Returns O(t+dt). Set heisenberg=true for evolving an observable in the heisenberg picture. If heisenberg=false then it is assumed that O is a density matrix.\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)","category":"page"},{"location":"documentation/#PauliStrings.rk4-Tuple{Function, Operator, Real, Real}","page":"Documentation","title":"PauliStrings.rk4","text":"rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)\n\nSingle step of Runge–Kutta-4 with time dependant Hamiltonian. H is a function that takes a number (time) and returns an operator.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Tools","page":"Documentation","title":"Tools","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"compress(o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.compress-Tuple{Operator}","page":"Documentation","title":"PauliStrings.compress","text":"compress(o::Operator)\n\nAccumulate repeated terms and remove terms with a coeficient smaller than 1e-20\n\n\n\n\n\n","category":"method"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"op_to_strings(o::Operator)","category":"page"},{"location":"documentation/#PauliStrings.op_to_strings-Tuple{Operator}","page":"Documentation","title":"PauliStrings.op_to_strings","text":"op_to_strings(o::Operator)\n\ntakes an operator, return (coefs, strings) where coefs is a list of numbers and strings is a list of pauli string coefs[i] multiply strings[i]\n\n\n\n\n\n","category":"method"},{"location":"documentation/#Index","page":"Documentation","title":"Index","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"","category":"page"},{"location":"docstrings/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"docstrings/","page":"Index","title":"Index","text":"","category":"page"}]
}
