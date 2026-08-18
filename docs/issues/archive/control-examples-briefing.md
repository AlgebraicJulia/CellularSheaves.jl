# Control Example Briefing

## Purpose

This briefing collects candidate control examples for future CellularSheaves.jl documentation and feature work around:

```math
\dot{x}(t) = A x(t) + B u(t),
```

zero-order-hold discretization, trajectory sheaves, harmonic-extension-based feasible trajectory spaces, and convex quadratic optimal control.

Selection criteria:

1. standard benchmark in broader mechanical engineering / robotics / vehicle control literature
2. intuitive story for documentation
3. compatible with linear or standard linearized control models
4. naturally extensible to networked / consensus / fleet / formation versions later

## Recommended shortlist

### 1. Double integrator / point mass

**Brief description.** A point mass moving on a line under an applied force or acceleration. This is the canonical second-order mechanical control model and the cleanest first LQR example.

**Typical variables.**

```math
x = \begin{bmatrix} p \\ \dot{p} \end{bmatrix},
\qquad
u = F
\quad\text{or}\quad
u = a,
```

```math
A = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix},
\qquad
B = \begin{bmatrix} 0 \\ 1 \end{bmatrix}.
```

**Why it is a good fit now.** Minimal state dimension, closed-form intuition, easy phase-plane plots, and perfect for validating trajectory-space and quadratic-cost code.

**Why it is a good fit later for consensus.** Vehicle platoons, rendezvous, and formation-control models often reduce to networks of double integrators.

**References.**

- Russ Tedrake, *Underactuated Robotics*, LQR chapter, MIT.
- Reza Olfati-Saber, J. Alex Fax, Richard M. Murray, “Consensus and Cooperation in Networked Multi-Agent Systems,” *Proceedings of the IEEE*, 2007.

### 2. Vehicle platoon / longitudinal car-following

**Brief description.** A convoy of vehicles moving in one dimension, with each vehicle modeled by position and velocity and controlled through throttle/brake or commanded acceleration.

**Typical variables.**

```math
x_i = \begin{bmatrix} p_i \\ \dot{p}_i \end{bmatrix},
\qquad
u_i = a_i.
```

**Why it is a good fit now.** It is still mathematically close to the double integrator, but already tells a concrete “machines and vehicles” story.

**Why it is a good fit later for consensus.** This is one of the most direct consensus extensions: spacing regulation, string stability, and distributed coordination on a path graph.

**References.**

- Devendra Swaroop, J. Karl Hedrick, “String Stability of Interconnected Systems,” *IEEE Transactions on Automatic Control*, 1996.
- Yang Zheng, Shengbo Eben Li, Keqiang Li, Fei-Yue Wang, “Stability and Scalability of Homogeneous Vehicular Platoon,” *IEEE Transactions on Intelligent Transportation Systems*, 2016.

### 3. Cart-pole / inverted pendulum on a cart

**Brief description.** A cart moving on a line with an actuated horizontal force, balancing an inverted pendulum near the upright equilibrium.

**Typical variables.**

```math
x =
\begin{bmatrix}
q \\
\theta \\
\dot{q} \\
\dot{\theta}
\end{bmatrix},
\qquad
u = F.
```

**Why it is a good fit now.** Famous unstable benchmark, visually intuitive, and a great example of linearization plus LQR on a real mechanical system.

**Why it is a good fit later for consensus.** Less direct than platoons, but arrays of balancing agents or coupled underactuated units give a natural networked extension.

**References.**

- Russ Tedrake, *Underactuated Robotics*, cart-pole chapters.
- Mark W. Spong, Seth Hutchinson, M. Vidyasagar, *Robot Modeling and Control*, Wiley, 2006.

### 4. Planar quadrotor / hover linearization

**Brief description.** A rigid body moving in a vertical plane with two thrust inputs. Around hover, the system linearizes to a standard controllable LTI model.

**Typical variables.**

```math
x =
\begin{bmatrix}
y \\
z \\
\phi \\
\dot{y} \\
\dot{z} \\
\dot{\phi}
\end{bmatrix},
\qquad
u =
\begin{bmatrix}
u_1 \\
u_2
\end{bmatrix}.
```

**Why it is a good fit now.** Strong robotics flavor, good plots, and a clear continuous-time-to-discrete-time story.

**Why it is a good fit later for consensus.** Quadrotor swarms and formation flight are standard multi-agent benchmarks.

**References.**

- Russ Tedrake, *Underactuated Robotics*, quadrotor examples.
- Robert Mahony, Vijay Kumar, Peter Corke, “Multirotor Aerial Vehicles: Modeling, Estimation, and Control of Quadrotor,” *IEEE Robotics & Automation Magazine*, 2012.

### 5. Lateral bicycle model / lane-keeping vehicle

**Brief description.** A vehicle lateral-dynamics model for steering control, usually linearized at constant forward speed and used for lane keeping or path tracking.

**Typical variables.**

```math
x =
\begin{bmatrix}
e \\
\dot{e} \\
\psi_e \\
\dot{\psi}_e
\end{bmatrix},
\qquad
u = \delta_f.
```

**Why it is a good fit now.** Real vehicle-control benchmark with a compelling engineering story and a standard LTI approximation.

**Why it is a good fit later for consensus.** Fleets of autonomous cars following shared lanes or formation trajectories make a natural networked extension.

**References.**

- Rajesh Rajamani, *Vehicle Dynamics and Control*, Springer, 2012.
- Jason Kong, Mark Pfeiffer, Georg Schildbach, Francesco Borrelli, “Kinematic and Dynamic Vehicle Models for Autonomous Driving Control Design,” IEEE IV workshop notes, 2015.

## Additional strong candidates

### 6. Mass-spring-damper

**Brief description.** A mass attached to a spring and damper, driven by an external force.

**Typical variables.**

```math
x =
\begin{bmatrix}
q \\
\dot{q}
\end{bmatrix},
\qquad
u = F.
```

**Why it is a good fit now.** Clean mechanical interpretation, standard second-order benchmark, and natural energy story.

**Why it is a good fit later for consensus.** Chains or networks of masses and springs map naturally to graph-coupled systems.

**References.**

- B. D. O. Anderson, John B. Moore, *Optimal Control: Linear Quadratic Methods*, 1990.
- Hassan K. Khalil, *Nonlinear Systems*, example sections on linear second-order systems.

### 7. DC motor / servo drive

**Brief description.** A motor or servo with electrical input and rotational mechanical output, used all over machine control.

**Typical variables.**

```math
x =
\begin{bmatrix}
\theta \\
\omega
\end{bmatrix}
\quad\text{or}\quad
\begin{bmatrix}
\theta \\
\omega \\
i
\end{bmatrix},
\qquad
u = V_a.
```

**Why it is a good fit now.** Strong “machine” example with direct industrial relevance and well-known LQR formulations.

**Why it is a good fit later for consensus.** Multi-axis servo coordination and robot-joint synchronization are natural follow-ons.

**References.**

- Karl J. Åström, Björn Wittenmark, *Computer-Controlled Systems*, 1997.
- Gene F. Franklin, J. David Powell, Abbas Emami-Naeini, *Feedback Control of Dynamic Systems*.

### 8. Satellite attitude / formation

**Brief description.** A rigid spacecraft or linearized relative-orbit model controlled by torques or small thrusts.

**Typical variables.**

```math
x =
\begin{bmatrix}
\phi \\
\theta \\
\psi \\
\dot{\phi} \\
\dot{\theta} \\
\dot{\psi}
\end{bmatrix},
\qquad
u =
\begin{bmatrix}
\tau_x \\
\tau_y \\
\tau_z
\end{bmatrix}.
```

**Why it is a good fit now.** Clean multi-input LTI example with geometric meaning and a long optimal-control literature.

**Why it is a good fit later for consensus.** Satellite formation flying is a standard coordinated-control application.

**References.**

- Peter C. Hughes, *Spacecraft Attitude Dynamics*, Dover.
- Kyle T. Alfriend et al., *Spacecraft Formation Flying*, Elsevier.

## Vetting notes

### Best immediate examples

1. **Double integrator**
   - best first documentation example
   - easiest to verify analytically
   - strongest foundation for later consensus variants

2. **Vehicle platoon**
   - best bridge from single-agent to networked / consensus control
   - directly relevant to future graph-based demos
   - still simple enough to explain briefly

3. **Cart-pole**
   - best unstable mechanical benchmark
   - classic control example everyone recognizes
   - good for showing linearization before consensus work

4. **Planar quadrotor**
   - best robotics-flavored example
   - good visuals and strong swarm extension story

5. **Mass-spring-damper**
   - best simple “machine / vibration” example
   - very clean consensus and coupled-network extension

### Probably weaker as first examples

- **DC motor / servo drive**: strong engineering relevance, but less visually compelling.
- **Satellite attitude**: elegant, but probably less immediate than robots or vehicles for the current docs.
- **Lateral bicycle model**: strong if the repo wants an autonomous-driving emphasis, but less fundamental than the double integrator or cart-pole.

## Suggested sequence for documentation

### Option A: simplest progression

1. single integrator
2. double integrator
3. vehicle platoon

This gives the clearest path into future consensus examples.

### Option B: broader mechanics progression

1. double integrator
2. cart-pole
3. mass-spring-damper chain

This gives a more classical mechanics flavor.

### Option C: robotics / vehicles progression

1. double integrator
2. planar quadrotor
3. vehicle platoon

This gives the strongest robotics-to-consensus story.

## My recommendation

If the next step is consensus-oriented, the best set to carry forward is:

1. **Double integrator**
2. **Vehicle platoon**
3. **Planar quadrotor** or **cart-pole**
4. **Mass-spring-damper chain**

That mix gives:

- one canonical analytic benchmark
- one explicitly network-ready vehicle example
- one visually compelling robotics or mechanics example
- one mechanically coupled system that already looks graph-like
