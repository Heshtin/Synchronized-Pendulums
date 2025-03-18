# Weak/Measure Synchronization in Coupled Pendulum Systems Suspended From a Common Beam

**Synchronization** is a ubiquitous form of cooperative behavior found in nature—ranging from the synchronized oscillations of metronomes on a common platform, to the integer multiple resonances in Jupiter’s moons, to the collective flashing of fireflies, and even to coordinated biochemical processes in organisms. While the underlying mechanisms can be intricate—for example, the minute feedback adjustments in fireflies or the biochemical pathways in living organisms—the core principle behind simpler synchronized systems (e.g., metronomes or orbital resonances of moons) is often surprisingly straightforward. These simpler dynamical systems have only a few degrees of freedom, making it easier to analyze and observe clear synchronized behavior.

Christiaan Huygens [1], in 1665, was the first to document such natural synchronization phenomena when he observed that two pendulum clocks mounted on a common beam always settled into oscillations at the same frequency and maintained a near-constant phase difference of \(\pi\). Later investigations showed that both **in-phase** and **anti-phase** synchronization can emerge in Huygens’ setup. For instance, Dilão [2] introduced a mechanism based on limit cycles in phase space that explains exact in-phase or anti-phase synchronization, while Czolczynski and co-workers [3] demonstrated experimentally and numerically that both in-phase and anti-phase synchronization can arise depending on initial conditions and system parameters. This strong or **complete synchronization** has been extensively studied under the name of Huygens’ synchronization.

A less-studied but more intriguing phenomenon is a **weaker form of synchronization**—also referred to as **measure synchronization**—in which two systems are not strictly phase-locked, yet their trajectories fill out the same region of phase space and share the same invariant measure. Numerically, this kind of synchronization can be observed in both regular and chaotic dynamical systems [4]. In principle, one can derive it from the action–angle formulation for regular systems, suggesting it is a broadly occurring manifestation of coherent evolution in Hamiltonian dynamics.

In this paper, we further analyze the phenomenon of measure synchronization as it appears within Huygens’ original system. Specifically, we consider two identical, non-dissipative pendulums, both of which are suspended from a common beam that is also free to move along its long axis in a frictionless (non-dissipative) manner. To exhibit measure synchronization, the two pendulums must be identical in mass. Only then does the natural symmetry of the system lead to orbits that occupy the same region of phase space, revealing this richer, weaker form of synchronization.

> **Note**: This project was inspired by a video by Veritasium:  
> [Veritasium: Coupled Pendulums Video](https://www.youtube.com/watch?v=t-_VPRCtiUg&t=85s)

**Weak synchronization,** also referred to as measure synchronization, is fascinating in coupled Hamiltonian systems because, under this form of synchronization, each subsystem occupies an identical region of phase space with the same invariant measure (i.e., the generalized volume in phase space). Even in our simple Huygens setup, this behavior can be surprisingly rich and beautiful due to the limited number of degrees of freedom, which makes it easier to observe and analyze.

---

## Huygens’ Synchronization Setup

A schematic of the Huygens synchronization setup can be found in various references (e.g., [link to a related study](https://royalsocietypublishing.org/doi/10.1098/rsos.170777)). In essence, two pendulums are suspended from a beam that can move along its horizontal axis. The equations of motion governing this system are:

A schematic of the Huygens synchronization setup can be found in various references (e.g., [link to a related study](https://royalsocietypublishing.org/doi/10.1098/rsos.170777)). In essence, two pendulums are suspended from a beam that can move along its horizontal axis. The equations of motion governing this system are:

A schematic of the Huygens synchronization setup can be found in various references (e.g., [link to a related study](https://royalsocietypublishing.org/doi/10.1098/rsos.170777)). In essence, two pendulums are suspended from a beam that can move along its horizontal axis. The equations of motion governing this system are:

$$
\begin{aligned}
m{l_1}^2\,\ddot{\theta}_1
&+ m g l_1 \sin(\theta_1)
+ m l_1 \,\ddot{x}\,\cos(\theta_1)
= 0,\\
m{l_2}^2\,\ddot{\theta}_2
&+ m g l_2 \sin(\theta_2)
+ m l_2 \,\ddot{x}\,\cos(\theta_2)
= 0,\\
(M + 2m)\,\ddot{x}
&+ m l_1 \Bigl(\ddot{\theta}_1 \cos(\theta_1) - \dot{\theta}_1^2 \sin(\theta_1)\Bigr)
+ m l_2 \Bigl(\ddot{\theta}_2 \cos(\theta_2) - \dot{\theta}_2^2 \sin(\theta_2)\Bigr)
= 0.
\end{aligned}
$$

Now, with a linear drag force with drag coefficient \(C_x\) and a spring with spring constant \(K\) attached to the rod, only the equation of motion of \(x\) changes:

$$
(M+2m)\,\ddot{x}
+ m l_1 \Bigl(\ddot{\theta}_1 \cos(\theta_1) - \dot{\theta}_1^2 \sin(\theta_1)\Bigr)
+ m l_2 \Bigl(\ddot{\theta}_2 \cos(\theta_2) - \dot{\theta}_2^2 \sin(\theta_2)\Bigr)
+ C_x\,\dot{x}
+ K\,x
= 0.
$$

Combining the first two equations of motion for \(\ddot{\theta}_i\) with the one for \(x\), we get:

$$
\ddot{x}
= \frac{
\frac{g}{2}\,\sin\bigl(2\theta_1\bigr)
+ \frac{g}{2}\,\sin\bigl(2\theta_2\bigr)
+ l_1\,\dot{\theta}_1^2\,\sin(\theta_1)
+ l_2\,\dot{\theta}_2^2\,\sin(\theta_2)
- C_x\,\dot{x}
- K\,x
}{
\frac{M}{m}
+ 2
- \cos^2(\theta_1)
- \cos^2(\theta_2)
}.
$$


---

## Numerical Simulations

To explore this system, we employ a fourth-order Runge–Kutta (RK4) integration scheme for the above set of ordinary differential equations (ODEs). By varying initial conditions, we can observe how the system transitions from the well-known strong (in-phase or anti-phase) synchronization to the more subtle measure synchronization. 

---

## References

1. C. Huygens, *Horologium Oscillatorium* (Maguy, Paris, 1673).  
2. R. Dilão, “Huygens synchronization of two clocks,” *Chaos* **19**, 023118 (2009).  
3. K. Czolczynski, P. Perlikowski, A. Stefanski, and T. Kapitaniak, “Huygens’ synchronization,” *International Journal of Bifurcation and Chaos* **21**, 2047 (2011).  
4. [Various references on measure synchronization in Hamiltonian systems]
