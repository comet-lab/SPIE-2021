# Hale-Bopp

What this code does
-------------------
This Matlab program simulates the exploration of the vocal folds with an endoscope and steerable wrist combination. It is based on code that simulated exploration of the middle ear with a trans-Eustachian endoscope.
The algorithm works as follows: first, it runs a sampling-based path planning algorithm (RRT) to generate a set of reachable points within the larynx. It then casts rays from each of these points to generate a map of the visible surface.

For further information see *Chiluisa et al. ISMR 2020*.

UPDATE LINK IF MOVING IT TO PERSONAL ACCOUNT CHANGES URL

For documentation on how to use this repository: https://wpi-comet-lab.gitbook.io/hale-bopp-documentation/-MC7XkSzIvo1U7f1XOjm/
 
How to run this code (also in the Gitbook)
--------------------
  **Step 1: Run a Simulation**  
  From inside Matlab, open the `src` directory and run `runSimulation.m`:
  
  If successful, three things will appear:
  1. In the command window, text detailing the results of RRT and ray-casting
  2. An animation figure with a model of the endoscope in the larynx model on the left & scatterplots on the right
  3. A figure with copies of the same model displaying a visibility map
  
  **Step 2: Have fun!**
