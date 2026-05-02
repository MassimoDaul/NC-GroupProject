# **SPH Group Project**

### Massimo Daul & Caleb Ashes

We built the initial simulation using Julia and an llm to help us debug and work faster. We built the visualizer in Julia as well (re the assignment instructions). We upload the simulation file, the visualizer file, and an .mp4 that records our simulation with smooth initial conditions. 

Our responses to the reflection questions are here:

1) Changing Initial Parameters

We tried changing a few of the main parameters:

-  increasing the time step (dt) made the simulation unstable really quickly. The particles started moving too fast, jittering, and sometimes flying out of the box. This makes sense because we’re taking bigger steps in time, so the integration becomes less accurate (not as smooth)

- increasing the stiffness (k) made the fluid feel more “solid” since pressure forces got stronger, but it also made things less stable. Bigger forces mean the system is harder to control numerically.

- increasing the viscosity (mu) made everything smoother. The particles didn’t move as chaotically and the fluid looked thicker, kind of like syrup.

- increasing the smoothing length (h / radius for neighborhood search) made the simulation smoother overall because each particle interacts with more neighbors. But the downside is that the fluid loses some detail and looks more blurry.

   
