# Overview
The goal of the project is to find the kidnapped vehicle which has been transported to a new location. To achieve this, a 2 dimensional particle filter is being implemented in c++ which uses a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data as inputs.

## Particle Filter Algorithm

The flowchart below represents the steps of the particle filter algorithm as well as its inputs.

![](images/flowchart.png)

1. Initialize the car position for all particles. If already initialized, predict the vehicle's next state plus all the particles from the previous data and estimated speed and yaw rate.
2. Get the landmark observations data from the simulator.
3. Update the particle weights and resample particles.
4. Calculate and output the average weighted error of the particle filter over all time steps so far.

## Inputs to the Particle Filter
The inputs to the particle filter are in the `data` directory. 

### The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

### All other data the simulator provides, such as observations and controls.

> * Map data provided by 3D Mapping Solutions GmbH.

## Success Criteria
If your particle filter passes the current grading code in the simulator (you can make sure you have the current version at any time by doing a `git pull`), then you should pass! 

The things the grading code is looking for are:


1. **Accuracy**: your particle filter should localize vehicle position and yaw to within the values specified in the parameters `max_translation_error` and `max_yaw_error` in `src/main.cpp`.

2. **Performance**: your particle filter should complete execution within the time of 100 seconds.