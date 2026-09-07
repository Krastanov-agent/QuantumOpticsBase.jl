using QuantumOptics
using CairoMakie


b = SpinBasis(1//2)
psi = (spinup(b) + spindown(b)) / sqrt(2)
fig, ax, plot = blochsphereplot(psi;
    axis=(title="Pure spin along x", xlabel="⟨σx⟩", ylabel="⟨σy⟩", zlabel="⟨σz⟩"),
    figure=(size=(640, 480),))
save("bloch-pure.png", fig) #hide
nothing #hide


b = SpinBasis(1//2)
rho = 0.7dm(spinup(b)) + 0.3dm(spindown(b))
fig = Figure(size=(640, 480))
ax = Axis3(fig[1, 1]; title="70% spin up, 30% spin down",
    xlabel="⟨σx⟩", ylabel="⟨σy⟩", zlabel="⟨σz⟩", aspect=:equal)
blochsphereplot!(ax, rho)
save("bloch-mixed.png", fig) #hide
nothing #hide


fig = with_theme(theme_dark(); BlochSpherePlot=(spherecolor=(:gray65, 0.2),
        wireframecolor=(:gray80, 0.6), wireframewidth=1.5,
        sphereresolution=(16, 8), color=:white)) do
    b = SpinBasis(1//2)
    psi = (spinup(b) + im * spindown(b)) / sqrt(2)
    fig, ax, plot = blochsphereplot(psi;
        axis=(title="Pure spins along y and z", xlabel="⟨σx⟩", ylabel="⟨σy⟩", zlabel="⟨σz⟩"),
        figure=(size=(640, 480),))
    blochsphereplot!(ax, spinup(b); spherevisible=false, color=:orange)
    fig
end
save("bloch-dark.png", fig) #hide
nothing #hide


b = FockBasis(20)
fig, ax, plot = fockdistributionplot(coherentstate(b, 2);
    axis=(title="Coherent state, α = 2", xlabel="Occupation n", ylabel="P(n)"),
    figure=(size=(640, 400),))
save("fock-coherent.png", fig) #hide
nothing #hide


b = FockBasis(12, 4)
fig, ax, plot = fockdistributionplot(fockstate(b, 7);
    axis=(title="Number state in a basis starting at n = 4",
        xlabel="Occupation n", ylabel="P(n)", xticks=4:12),
    figure=(size=(640, 400),))
save("fock-offset.png", fig) #hide
nothing #hide


b = FockBasis(30)
psi = coherentstate(b, 2)
rho = 0.5dm(coherentstate(b, 1)) + 0.5dm(coherentstate(b, 3))
fig = Figure(size=(640, 400))
ax = Axis(fig[1, 1]; title="Coherent state and a classical mixture",
    xlabel="Occupation n", ylabel="P(n)", limits=(-0.5, 20.5, nothing, nothing))
fockdistributionplot!(ax, psi; dodge=1, n_dodge=2, label="α = 2")
fockdistributionplot!(ax, rho; dodge=2, n_dodge=2, label="Mixture of α = 1 and α = 3")
axislegend(ax)
save("fock-comparison.png", fig) #hide
nothing #hide


b = FockBasis(30)
x = range(-4, 5; length=151)
p = range(-4, 4; length=151)
fig, ax, plot = wignerplot(coherentstate(b, 1 + 0.5im), x, p;
    axis=(title="Coherent-state Wigner function", xlabel="Position x", ylabel="Momentum p", aspect=DataAspect()),
    figure=(size=(640, 480),))
Colorbar(fig[1, 2], plot; label="W(x, p)")
save("wigner-coherent.png", fig) #hide
nothing #hide


fig = with_theme(Theme(WignerPlot=(colorrange=(-1/pi, 1/pi),))) do
    b = FockBasis(10)
    x = p = range(-4, 4; length=151)
    fig, ax, plot = wignerplot(fockstate(b, 1), x, p;
        axis=(title="Single-photon Wigner function", xlabel="Position x", ylabel="Momentum p", aspect=DataAspect()),
        figure=(size=(640, 480),))
    Colorbar(fig[1, 2], plot; label="W(x, p)")
    fig
end
save("wigner-negative.png", fig) #hide
nothing #hide


b = FockBasis(30)
psi = normalize(coherentstate(b, 2) + coherentstate(b, -2))
x = range(-5, 5; length=201)
p = range(-3, 3; length=151)
fig = Figure(size=(640, 440))
ax = Axis(fig[1, 1]; title="Even cat state", xlabel="Position x", ylabel="Momentum p", aspect=DataAspect())
plot = wignerplot!(ax, dm(psi), x, p; colorrange=(-1/pi, 1/pi))
Colorbar(fig[1, 2], plot; label="W(x, p)")
save("wigner-cat.png", fig) #hide
nothing #hide


b = PositionBasis(-6, 6, 256)
psi = gaussianstate(b, 1, 2, 1)
fig, ax, plot = wavefunctionplot(psi;
    axis=(title="Gaussian position density", xlabel="Position x", ylabel="|ψ(x)|²"),
    figure=(size=(640, 400),))
save("wavefunction-position.png", fig) #hide
nothing #hide


b = PositionBasis(-6, 6, 256)
psi = gaussianstate(b, 0, 3, 1.5)
fig = Figure(size=(640, 400))
ax = Axis(fig[1, 1]; title="Complex Gaussian wavefunction",
    xlabel="Position x", ylabel="Wavefunction amplitude")
wavefunctionplot!(ax, psi; component=real, label="Re ψ(x)")
wavefunctionplot!(ax, psi; component=imag, label="Im ψ(x)")
axislegend(ax)
save("wavefunction-components.png", fig) #hide
nothing #hide


bx = PositionBasis(-32, 32, 512)
bp = MomentumBasis(bx)
psi = transform(bp, bx) * gaussianstate(bx, 0, 2, 1)
fig, ax, plot = wavefunctionplot(psi;
    axis=(title="Gaussian momentum density", xlabel="Momentum p", ylabel="|ψ(p)|²", limits=(-3, 7, nothing, nothing)),
    figure=(size=(640, 400),))
save("wavefunction-momentum.png", fig) #hide
nothing #hide
