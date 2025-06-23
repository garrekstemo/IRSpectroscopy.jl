using LVM
using GLMakie

mol = LVM.readlvm("data/sig_250619_172532.lvm")
cav = LVM.readlvm("data/sig_250619_185240.lvm")

ΔA_mol = -log10.(mol.on ./ mol.off)
ΔA_cav = -log10.(cav.on ./ cav.off)
diff_mol = mol.on .- mol.off
diff_cav = cav.on .- cav.off

fig = Figure(size = (500, 700))
DataInspector()

offset = 2
ax1 = Axis(fig[1, 1])
lines!(mol.wavenumber, ΔA_mol ./ maximum(ΔA_mol))
# lines!(mol.wavenumber, diff_mol ./ maximum(diff_mol))

lines!(mol.wavenumber, mol.off ./ maximum(-mol.off) .+ offset, label = "off")
lines!(mol.wavenumber, mol.on ./ maximum(-mol.on) .+ offset, label = "on")


ax2 = Axis(fig[2, 1])
lines!(cav.wavenumber, ΔA_cav ./ maximum(ΔA_cav))
# lines!(cav.wavenumber, diff_cav ./ maximum(diff_cav))

lines!(cav.wavenumber, cav.off ./ maximum(-cav.off) .+ offset, label = "off")
lines!(cav.wavenumber, cav.on ./ maximum(-cav.on) .+ offset, label = "on")


axislegend(ax1)
axislegend(ax2)

fig
