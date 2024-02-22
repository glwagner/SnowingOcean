using GLMakie
using GibbsSeaWater
using Printf

fig = Figure(size=(1200, 400))
axp = Axis(fig[1, 1], xlabel="Depth (m)", ylabel="Freezing temperature (ᵒC)")
axS = Axis(fig[1, 2], xlabel="Absolute salinity (g kg⁻¹)", ylabel="Freezing temperature (ᵒC)")
axa = Axis(fig[1, 3], xlabel="Absolute salinity (g kg⁻¹)", ylabel="Thermal expansion \n at freezing temperature (kg m⁻³ ᵒC⁻¹")

p₀ = 10.1325
p = p₀ .+ 0.0:10.0:1000

for Sᴬ = (0, 20, 30, 40)
    T = GibbsSeaWater.gsw_ct_freezing.(Sᴬ, p, 0)
    lbl = @sprintf("Sᴬ = %d g/kg", Sᴬ)
    lines!(axp, p, T, label=lbl)
end
    
Sᴬ = 0.0:0.1:50.0

for p = (0, 100, 1000)
    T = GibbsSeaWater.gsw_ct_freezing.(Sᴬ, p + p₀, 0)
    lbl = @sprintf("z = -%d m", p)
    lines!(axS, Sᴬ, T, label=lbl)
end

for p = (0, 100, 1000)
    T★ = GibbsSeaWater.gsw_ct_freezing.(Sᴬ, p + p₀, 0)
    α = gsw_alpha.(Sᴬ, T★, p + p₀)
    lbl = @sprintf("z = -%d m", p)
    lines!(axa, Sᴬ, α, label=lbl)
end
 
axislegend(axp, position=:rt)
axislegend(axS, position=:lb)
axislegend(axa, position=:rb)

display(fig)

