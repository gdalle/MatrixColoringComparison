### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ da55551a-228e-11ef-3448-ef1ce3658ba7
begin
	using Pkg
    Pkg.activate(mktempdir())
	Pkg.add([
	    Pkg.PackageSpec("CairoMakie"),
		Pkg.PackageSpec("Chairmarks"),
		Pkg.PackageSpec(url="https://github.com/exanauts/ColPack.jl", rev="gd/super_upgrade_colpack"),
		Pkg.PackageSpec("DataFrames"),
		Pkg.PackageSpec("DataFramesMeta"),
		Pkg.PackageSpec("LinearAlgebra"),
		Pkg.PackageSpec("PlutoUI"),
		Pkg.PackageSpec("ProgressLogging"),
		Pkg.PackageSpec("SparseMatrixColorings"),
		Pkg.PackageSpec("SparseArrays"),
		Pkg.PackageSpec("StableRNGs"),
		Pkg.PackageSpec("Statistics"),
		Pkg.PackageSpec("Test"),
    ])
	
    using CairoMakie
	using Chairmarks
	using ColPack
	using DataFrames
	using DataFramesMeta
	using LinearAlgebra
	using PlutoUI
	using ProgressLogging
	using SparseMatrixColorings
	using SparseMatrixColorings: NaturalOrder
	using SparseArrays
	using StableRNGs
	using Statistics
	using Test
end

# ╔═╡ a36b9ca9-a9e4-44df-8eca-991c2f3aa2e0
md"""
# Matrix coloring comparison
"""

# ╔═╡ 6e1fed93-3b43-41eb-9993-ac8ced7244e7
md"""
A comparison between matrix coloring packages in Julia.

At the moment we include:

- [SparseMatrixColorings.jl](https://github.com/gdalle/SparseMatrixColorings.jl)
- [ColPack.jl](https://github.com/exanauts/ColPack.jl)
"""

# ╔═╡ 1770b4f7-4021-4345-b191-2f92ce41deed
TableOfContents()

# ╔═╡ 867c780f-beed-4b10-a5dc-eb1385ef3941
md"""
## Setup
"""

# ╔═╡ 54524afe-643d-4c54-bb3c-95de77c2c7bf
md"""
### SparseMatrixColorings
"""

# ╔═╡ 733a96e6-e5c9-435b-861b-725ce0a5a48d
function smc_column_coloring(J)
	return column_coloring(J, GreedyColoringAlgorithm(NaturalOrder()))
end

# ╔═╡ 0e776e0b-23c9-484f-8cd5-457c9f9dc088
function smc_row_coloring(J)
	return row_coloring(J, GreedyColoringAlgorithm(NaturalOrder()))
end

# ╔═╡ 434a9de0-45d5-4c42-8e23-5ea7ab0437a1
function smc_symmetric_coloring(H)
	return symmetric_coloring(H, GreedyColoringAlgorithm(NaturalOrder()))
end

# ╔═╡ 1399b39f-1311-47c8-b474-b33a8955017f
md"""
### ColPack
"""

# ╔═╡ 75e1f722-dfdd-486d-8fd1-bb62909cd6e5
function colpack_column_coloring(J)
	method = "COLUMN_PARTIAL_DISTANCE_TWO"
	order = "NATURAL"
    coloring = ColPackPartialColoring(J, method, order)
    return get_colors(coloring)
end

# ╔═╡ 0f1bbfde-dede-4a8a-8914-9046a4673b88
function colpack_row_coloring(J)
    method = "ROW_PARTIAL_DISTANCE_TWO"
	order = "NATURAL"
    coloring = ColPackPartialColoring(J, method, order)
    return get_colors(coloring)
end

# ╔═╡ 18375b3b-a8b9-4d6d-8605-10525e5e2614
function colpack_symmetric_coloring(H)
    method = "STAR"
	order = "NATURAL"
    coloring = ColPackColoring(H, method, order)
    return get_colors(coloring)
end

# ╔═╡ 185bba7c-bde5-4787-a99f-7a6cce6130a1
md"""
### Misc
"""

# ╔═╡ 4f6561fc-15e0-4a28-9ae2-9c4979244222
d_values_small(n) = filter(<=(n), [5,10,20])

# ╔═╡ 7d6c9868-b770-4b2e-a4c2-fd55b2a61bf5
d_values(n) = filter(<=(n), [5,10,20,50,100])

# ╔═╡ b92ce3d8-0a67-4433-b722-031156b52cac
linestyles = [:solid, :dash, :dashdot, :dashdotdot, :dot]

# ╔═╡ 484036b2-8aba-476a-9336-2ddab9d1f23a
markers = [:circle, :star5, :rect, :utriangle, :cross]

# ╔═╡ ccace506-dcb4-4fa4-bc02-280c0c25f60b
md"""
## Correctness
"""

# ╔═╡ c9eae8aa-2163-4290-b479-91cf194be563
md"""
### Jacobian
"""

# ╔═╡ 86948d71-6933-422d-9345-d9e5f51f60ab
@testset "Column coloring" begin
	@testset "n=$n - d=$d" for n in 10 .^ (1:3), d in d_values_small(n)
		J = sprand(n, n + 1, d/n)
		color1 = smc_column_coloring(J)
		color2 = colpack_column_coloring(J)
		@test_broken color1 == color2
	end
end

# ╔═╡ 45220b76-353f-4c23-b2d9-aacff0e9a3c0
@testset "Row coloring" begin
	@testset "n=$n - d=$d" for n in 10 .^ (1:3), d in d_values_small(n)
		J = sprand(n, n + 1, d/n)
		color1 = smc_row_coloring(J)
		color2 = colpack_row_coloring(J)
		@test_broken color1 == color2
	end
end

# ╔═╡ 1854bc5c-70a9-4543-86e2-f81dc232959e
md"""
### Hessian
"""

# ╔═╡ 5b93fc63-435b-4d98-8b01-bea6ad2a84f5
@testset "Symmetric coloring" begin
	@testset "n=$n - d=$d" for n in 10 .^ (1:3), d in d_values_small(n)
		H = sparse(Symmetric(sprand(n, n, d/n)))
		color1 = smc_symmetric_coloring(H)
		color2 = colpack_symmetric_coloring(H)
		@test_broken color1 == color2
	end
end

# ╔═╡ b56392cc-3575-4c02-a672-c434263d803c
md"""
## Benchmarks
"""

# ╔═╡ 8e4a280d-9d38-447d-be1f-ed98ce856d83
function plot_comparison(data; matrix::Symbol)
	
	data = @orderby(data, :package, :matrix, :n, :d)
	smc_data = @rsubset(data, :package=="SparseMatrixColorings")
	cp_data = @rsubset(data, :package=="ColPack")
	ymax = 1.05 * maximum(cp_data[!, :time_min] ./ smc_data[!, :time_min])
	
	fig = Figure()
	ax = Axis(
		fig[1, 1], title="$matrix coloring",
		xlabel="dimension (number of rows & columns)",
		ylabel="speedup SMC / ColPack",
		xscale=log10,
		limits=(nothing, (0, ymax)),
	)

	hl = hlines!(ax, [1.0], linewidth=2, color=:black)
	for (i, d) in enumerate(sort(unique(data[!, :d])))
		smc = @rsubset(data, :d==d && :package=="SparseMatrixColorings")
		cp = @rsubset(data, :d==d && :package=="ColPack")
		sl = scatterlines!(
			smc[!, :n], cp[!, :time_min] ./ smc[!, :time_min];
			label="$d nz/col", linestyle=linestyles[i], marker=markers[i],
			markersize=12,
		)
	end
	band!(
		ax, [minimum(data[!, :n]), maximum(data[!, :n])], [0, 0], [1, 1],
		color=(:black, 0.1), label="ColPack\nfaster here"
	)
	Legend(fig[1, 2], ax)
	return fig
end

# ╔═╡ d1980626-0f1d-4af8-b7ba-1c798959b024
function benchmark_coloring(; n_values::AbstractVector{Int}, matrix::Symbol)
	data = DataFrame()
	nd_pairs = collect(
		(n, d)
		for n in sort(n_values; rev=true)
		for d in sort(d_values(n); rev=true)
	)
	@progress for (n, d) in nd_pairs
		p = d/n
		for package in ("SparseMatrixColorings", "ColPack")
			bench = if package == "SparseMatrixColorings" && matrix == :Jacobian
				@be sprand(n, n, p) smc_column_coloring seconds=5 evals=1 samples=100
			elseif package == "SparseMatrixColorings" && matrix == :Hessian
				@be sparse(Symmetric(sprand(n, n, p))) smc_symmetric_coloring seconds=5 evals=1 samples=100
			elseif package == "ColPack" && matrix == :Jacobian
				@be sprand(n, n, p) colpack_column_coloring seconds=5 evals=1 samples=100
			elseif package == "ColPack" && matrix == :Hessian
				@be sparse(Symmetric(sprand(n, n, p))) colpack_symmetric_coloring seconds=5 evals=1 samples=100
			end
			row = (;
				package=string(package),
				matrix=string(matrix),
				n=n, d=d,
				evals=Int(minimum(bench).evals),
				samples=length(bench.samples),
				time_min=minimum(bench).time,
				time_median=median(bench).time,
				time_q25=quantile(bench, 0.25).time,
				time_q75=quantile(bench, 0.75).time,
			)
			push!(data, row)
		end
	end
	return data
end

# ╔═╡ 6b5deada-4836-498b-aeb2-6817274733f9
md"""
### Jacobian
"""

# ╔═╡ 6f2a0f6e-068f-434e-af81-8a5b89b37e01
data_jacobian = benchmark_coloring(
	n_values = floor.(Int, 10 .^ (1:0.3:4)),
	matrix=:Jacobian,
)

# ╔═╡ 75824976-a2f9-4c5e-80b3-de5bae5ddaa3
plot_comparison(data_jacobian; matrix=:Jacobian)

# ╔═╡ dd62ca86-8a27-4628-8b66-8626ac7ae006
md"""
### Hessian
"""

# ╔═╡ 1011abc1-0877-462b-9cff-5dee8be097e0
data_hessian = benchmark_coloring(
	n_values=floor.(Int, 10 .^ (1:0.3:4)),
	matrix=:Hessian
)

# ╔═╡ 8155f0dd-1843-487b-96ff-a5c4cbc6db7c
plot_comparison(data_hessian; matrix=:Hessian)

# ╔═╡ Cell order:
# ╟─a36b9ca9-a9e4-44df-8eca-991c2f3aa2e0
# ╟─6e1fed93-3b43-41eb-9993-ac8ced7244e7
# ╠═da55551a-228e-11ef-3448-ef1ce3658ba7
# ╠═1770b4f7-4021-4345-b191-2f92ce41deed
# ╟─867c780f-beed-4b10-a5dc-eb1385ef3941
# ╟─54524afe-643d-4c54-bb3c-95de77c2c7bf
# ╠═733a96e6-e5c9-435b-861b-725ce0a5a48d
# ╠═0e776e0b-23c9-484f-8cd5-457c9f9dc088
# ╠═434a9de0-45d5-4c42-8e23-5ea7ab0437a1
# ╟─1399b39f-1311-47c8-b474-b33a8955017f
# ╠═75e1f722-dfdd-486d-8fd1-bb62909cd6e5
# ╠═0f1bbfde-dede-4a8a-8914-9046a4673b88
# ╠═18375b3b-a8b9-4d6d-8605-10525e5e2614
# ╟─185bba7c-bde5-4787-a99f-7a6cce6130a1
# ╠═4f6561fc-15e0-4a28-9ae2-9c4979244222
# ╠═7d6c9868-b770-4b2e-a4c2-fd55b2a61bf5
# ╠═b92ce3d8-0a67-4433-b722-031156b52cac
# ╠═484036b2-8aba-476a-9336-2ddab9d1f23a
# ╟─ccace506-dcb4-4fa4-bc02-280c0c25f60b
# ╟─c9eae8aa-2163-4290-b479-91cf194be563
# ╠═86948d71-6933-422d-9345-d9e5f51f60ab
# ╠═45220b76-353f-4c23-b2d9-aacff0e9a3c0
# ╟─1854bc5c-70a9-4543-86e2-f81dc232959e
# ╠═5b93fc63-435b-4d98-8b01-bea6ad2a84f5
# ╟─b56392cc-3575-4c02-a672-c434263d803c
# ╠═8e4a280d-9d38-447d-be1f-ed98ce856d83
# ╠═d1980626-0f1d-4af8-b7ba-1c798959b024
# ╟─6b5deada-4836-498b-aeb2-6817274733f9
# ╠═6f2a0f6e-068f-434e-af81-8a5b89b37e01
# ╠═75824976-a2f9-4c5e-80b3-de5bae5ddaa3
# ╟─dd62ca86-8a27-4628-8b66-8626ac7ae006
# ╠═1011abc1-0877-462b-9cff-5dee8be097e0
# ╠═8155f0dd-1843-487b-96ff-a5c4cbc6db7c
