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

The biggest difference I noticed was that the cubic spline kernel produced smoother particle motion. The fluid looked less noisy overall, especially in regions where particles were interacting heavily. The pressure forces also seemed slightly more stable, with fewer sudden spikes in movement.

However, the simulation also looked a bit more “blurred” compared to the original kernel. Fine details in the fluid motion were less sharp, which makes sense since the cubic spline spreads interactions more gradually across neighbors.

In terms of stability, the cubic spline kernel felt slightly more forgiving. The simulation didn’t break as quickly when parameters like dt or k were increased, suggesting that the smoother gradients help reduce numerical instability.

Overall, switching kernels changed the balance between detail and stability. The original kernels gave sharper structure but more noise, while the cubic spline produced smoother and more stable behavior at the cost of some detail.
