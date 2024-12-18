using CairoMakie
using Unitful
using StatsBase
using ColorSchemeTools
using ColorTypes

t = 30u"s"
a = final_ions(30u"s", 1000000, dNi_excitation, gaussian_sample, 10^8)
H = fit(Histogram, (ustrip.(getindex.(a,1)), ustrip.(getindex.(a,3))), nbins=5000)

scheme = make_colorscheme([RGB(0,0,0), RGB(0.5,0,1)], 500)

fig = Figure()
ax = Axis(fig[1,1],
    title = "Simulated ion cloud at t=$t",
    xlabel = "Injection Direction (km)",
    ylabel = "Magnetic Field Direction (km)",
    aspect = 1
)

heatmap!(ax, H, colormap=scheme)

fig