# Rocket trajectory with Gravity Turn
This is a simplified model for a rocket trajectory going into LEO (Low Earth Orbit) in python, taking as reference the Electron from RocketLab.

From the Figure 1, the trajectory in relation to the Earth, for simplicity as said, we will ignore the Earth´s rotation and write the equations respect to a non rotating earth.

![image](https://github.com/camilo-zuluaga/rocket-trajectory/assets/103959811/c3792283-d3c5-4509-9c73-dc513fac440e)

<div align="center">
  <strong>Figure: 1.</strong> Representation of a gravity turn
</div>

# Differential Equations

_T direction:_

$$\Large \dot{v} = -\frac{T}{m} - g \cos(\psi) + \frac{1}{2m} \rho A C_d v^2$$

_N direction:_

$$\Large \dot{\phi} = -\frac{g}{v} \sin(\psi)
$$
_From geometry:_
$$\Large \dot{h} = v \cos(\psi)
$$
$$\Large \dot{\theta} = \frac{v \sin(\psi)}{R_e + h}
$$
$$\Large \psi = \dot{\phi} - \dot{\theta}
$$

# Modeling enviroment
_Drag Force_

$$\Large D = \frac{1}{2} \rho v^2 C_D A$$

_Atmospheric Density_

$$\Large \rho = \rho_0 e^{-\frac{h}{H}}$$

_Gravity Acceleration_

$$\Large g = \frac{g_0}{\left(1 + \frac{h}{R_e}\right)^2}$$

# Vehicle Parameters

| LiftOff Mass (Kg) | First Stage (kN) | Second Stage (kN) | Payload (Kg) |
| :---------------: | :--------------: | :---------------: | :----------: |
|       13500       |     218.982      |       29.22       |     300      |

# Results

![Results](https://github.com/camilo-zuluaga/rocket-trajectory/assets/103959811/7e661874-e2e9-47ec-9a08-b6c353924b91)

# Disclaimer

This is an early version and some values havent been modified, next updates will make this in a more realistic orbit for the Electron Rocket.

