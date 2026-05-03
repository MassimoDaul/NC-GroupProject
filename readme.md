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

When we switched from the (more efficient) spatial hashing approach to the (less efficient) naïve neighbor search, the performance was a lot worse. This is clear just from the quadratic runtime, as opposed to the (quasi-)linear runtime of spatial hashing. 

**5) Changing Integration Scheme**

Finally, when we switched from the Euler integration scheme to the Euler-Cromer integration scheme, we noticed that short-time behavior was quite similar while long-time behavior was different. The Euler-Cromer scheme was more stable for longer time intervals, and we can easily explain this if we think of the scheme as a matrix equation. 

In the Euler scheme, the eigenvalues of the update matrix are bigger than 1 (in magnitude); in the Euler-Cromer scheme, the eigenvalues are exactly 1 (in magnitude), which explains why it is more stable for longer time intervals. 
