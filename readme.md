# **SPH Group Project**

### Massimo Daul & Caleb Ashes

We built the initial simulation using Julia and an llm to help us debug and work faster. We built the visualizer in Julia as well (re the assignment instructions). We upload the simulation file, the visualizer file, and an .mp4 that records our simulation with smooth initial conditions. 

Our responses to the reflection questions are here:

**1) Changing Initial Parameters**

We tried changing a few of the main parameters:

-  increasing the time step (dt) made the simulation unstable really quickly. The particles started moving too fast, jittering, and sometimes flying out of the box. This makes sense because we’re taking bigger steps in time, so the integration becomes less accurate (not as smooth)

- increasing the stiffness (k) made the fluid feel more “solid” since pressure forces got stronger, but it also made things less stable. Bigger forces mean the system is harder to control numerically.

- increasing the viscosity (mu) made everything smoother. The particles didn’t move as chaotically and the fluid looked thicker, kind of like syrup.

- increasing the smoothing length (h / radius for neighborhood search) made the simulation smoother overall because each particle interacts with more neighbors. But the downside is that the fluid loses some detail and looks more blurry.


**2) Using a Cubic Spline Kernal**

We tried replacing the original poly6 kernel combination with a cubic spline kernel and reran the simulation to compare the behavior. We used this:

```function W_cubic(r)
    q = r / h
    σ = 10 / (7 * pi * h^2)

    if 0 <= q < 1
        return σ * (1 - 1.5q^2 + 0.75q^3)
    elseif 1 <= q < 2
        return σ * 0.25 * (2 - q)^3
    else
        return 0.0
    end
end

function gradW_cubic(rvec)
    r = norm(rvec)
    if r == 0
        return zeros(2)
    end

    q = r / h
    σ = 10 / (7 * pi * h^2)

    if 0 <= q < 1
        dWdr = σ * (-3q + 2.25q^2) / h
    elseif 1 <= q < 2
        dWdr = -σ * 0.75 * (2 - q)^2 / h
    else
        dWdr = 0.0
    end

    return dWdr * (rvec / r)
end
```

We thought that the cubic spline kernel produced smoother particle motion. The fluid looked less noisy overall and there weren't really any spikes in the movement.  Simulation also looked more blurry compared to the original kernel. 


**3) Adding More Particles**

We changed the initial grid to use:
```
for x in dx:dx:0.6
    for y in dx:dx:0.8
```

With more particles, the fluid looked fuller and the motion had more detail. The dam break also covered more of the domain because there was simply more fluid at the beginning. It was cooler in general. The main downside was performance, which is a little obvious. It just took way longer to run. 

**4) Changing Neighborhood Search Algorithm**

The spatial hashing approach partitions the simulation domain into a grid of cells each of side length h (the smoothing length). For each particle, only the 3×3 block of surrounding cells needs to be examined to find all neighbors within the interaction radius h. Since the number of particles per cell is bounded by a constant k determined by the initial particle spacing, the total cost of the neighbor search is O(n · k) ≈ O(n), scaling linearly in the number of particles.

The naïve neighbor search instead checks every pair of particles (i, j) and tests whether their distance falls within h, yielding a cost of O(n²). If the number of particles is doubled, the naïve search requires four times as much computation whereas spatial hashing requires only twice as much. For the enlarged particle setup explored in TODO 8, this difference becomes severe — simulations completing in seconds under spatial hashing may take orders of magnitude longer under the naïve approach.

Switching to a naïve search has no effect on the physical accuracy of the simulation. Both methods identify exactly the same set of neighbors for each particle. The difference is purely computational: spatial hashing exploits the locality of SPH interactions, while the naïve search discards this structure at significant performance cost.

**5) Changing Integration Scheme**

This single difference has a decisive consequence for long-time energy behavior, which can be understood via eigenvalue analysis. For a simple harmonic oscillator, the exact solution traces a closed circle in phase space. Examining the eigenvalues of each scheme's one-step update matrix reveals the following. Under Euler the eigenvalues have modulus strictly greater than 1, so the numerical trajectory spirals outward in phase space and the total energy of the system grows monotonically over time. Under Euler-Cromer the eigenvalues have modulus exactly equal to 1, so the trajectory traces a slightly distorted but closed curve and energy oscillates around the true value with no long-time drift. This is the defining property of a symplectic integrator: it preserves the symplectic structure of Hamiltonian mechanics, guaranteeing bounded energy error over arbitrarily long simulations.

Applied to the SPH simulation, Euler causes particles to gradually accumulate kinetic energy over many time steps until pressure forces become unbalanced, eventually producing unphysical velocities and numerical blow-up. Short-time behavior appears similar between the two schemes since both are first-order accurate and the per-step energy error is small, but the instability compounds over time. Reducing dt slows the onset of instability but does not eliminate it, since the eigenvalue modulus remains greater than 1 for any nonzero dt under Euler. Euler-Cromer avoids this entirely and remains stable for the duration of the simulation provided dt is sufficiently small.
