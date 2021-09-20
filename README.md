# Physics_Collisions
### Classroom project realized with Isabella Sanna 
In a collider experiment (as LHC), particles collide in a region centered in the origin of a Cartesian coordinate system, where z axis is the direction of the beam.
The information coming from the vertex detector is later used to define where the interaction happened.

In this project we simulate the collisions and the particle transport in the space to reconstruct the vertex of interaction (z coordinate). By comparing the real and reconstructed positions of the vertex, we define the resolution of the detector under examination (around 180 micrometer).

Working hypothesis:
- detector made of two silicon pixel planes, coaxial with the beam axis;
- trajectory curvature due to magnetic field neglected;
- axial simmetry;
- vertex coordinates: z extracted from normal distribution with rms=5.3 cm, x and y extracted from norml distribution with rms=0.1 mm;
- Berilium beampipe: radius=3 cm, thickness=0.8 mm; first layer: radius=4 cm, thickness=0.2 mm; second layer: radius=7 cm;
- detector extention along z axis=27 cm (acceptance -1<&eta;<1 for particles produced in collisions with &sigma;<Zvertex<&sigma;);                    
- “fast” reconstruction with smearing of collision points: z -> 120 &mu;m r&phi; -> 30 &mu;m;
- multiple scattering in beampipe and layers.

By micol.olocco@gmail.com, isabella.sanna@edu.unito.it
