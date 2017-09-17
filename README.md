# **CarND: 2D Particle Filter**  [![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)
[//]: # (Image References)
[sample_gif]: ./img/pf_sample.gif
[process]: ./img/process_schema.png
[ctrv]: ./img/ctrv.png
[ctrv_yaw0]: ./img/ctrv_yaw0.png

The goal of this project is to implement a **2D particle filter** in C++ to localize a vehicle. Given the map, the particle filter is initialized with the help of a GPS measurement. The subsequent steps in the process will refine this estimate to localize our vehicle with the help of observation and control data.

![sample_gif]

The observations, control, and update state data are streamed with Udacity's [Self-Driving Car simulator](https://github.com/udacity/self-driving-car-sim).

## Process and Implementation

The flowchart bellow represents the steps and corresponding inputs to the particle filter algorithm:

![process]

### 0. The Particles

In this project, the vehicle is simplified with the **bycicle model**:
+ 2D motion
+ The front wheels of the car are connected to the rear wheels by a rigid beam with fixed length
+ The front two wheels are together, so they can be represented by one wheel. The same holds for the rear wheels.
+ The car is controlled with a steering angle of the simplified front wheel.

Bearing this in mind, the location of each particle will be determined by its `[x, y]` position and heading `theta`.

### 1. Initialization

The particles are instantiated with an initial **GPS estimate**. The particles are sampled from a normal distribution centered in the measured GPS coordinates, with variance equal to the GPS sensor noise.

This step is performed in the `ParticleFilter::init` method.

### 2. Prediction

Given the latest control inputs (yaw rate and velocity) of the vehicle, the position of each particle at the next time step is predicted by applying the **CTRV motion model** *(Constant Turn Rate and Velocity)*. For non-zero yaw rates, the update equations are:

![ctrv]

Wheareas if the yaw rate is zero:

![ctrv_yaw0]

This step is implemented in the `ParticleFilter::prediction` method.

*Note*: Gaussian noise should be added to the velocity and yaw rate inputs, in order to account for the uncertainty in the control inputs. However, as provided by Udacity, the method is expected to apply GPS sensor noise on `x`, `y` and `theta` for each particle. IMO this is erroneous, but I have yet decided to implementet it in order to fit the provided interface and meet the final grading criteria.

### 3. Data Association

For each particle, each one of the measured observations is associated to a map landmark via nearest neighbour search. In order to calculate observation-landmark distances, the observations are transformed from local (particle) to map (global) frame of reference.

```
for each particle :
    get landmarks in range
    for each observation:
        transform observation to map coordinates
        get nearest landmark in range to observation
```

This pseudocode snippet is implemented in `ParticleFilter::updateWeights`.

### 4. Update

The final weight of each particle is calculated as the product of each measurement's probability to correspond to its associated landmark. The probability factors are approximated by evaluating a Multivariate Gaussian distribution in `ParticleFilter::updateWeights`.

The mean of the distribution is the measurement's associated landmark position and the standard deviation is described by the initial uncertainty in the x and y ranges. The Multivariate-Gaussian is evaluated at the point of the transformed measurement's position. 
