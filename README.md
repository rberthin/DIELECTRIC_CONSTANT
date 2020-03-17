# DIELECTRIC_CONSTANT
Computation of the dielectric constant and error on this value. 
This script is in two parts -->

  1- The fortran script :
        Need epsilon.inpt (all the input parameters), data.inpt (MetalWalls input file, with the position of the step 0), 
        dipoles.out (MetalWalls ouput file), polarization.out (MetalWalls ouput file, total_dipole/volume).
        This script first rebuilt the box, compute total dipole (induced+permanent) on step 0 and then substracte this 
        value to step 0 in polarization.out. This difference is then substracted on each step. 
        The output file is total_polarization.dat that contain for each step sqrt(px**2+py**2+pz**2)
        
  2- The python script :
        Compute the correlation function of the polarization.
        Compute epsilon and error on epsilon.
        
The both codes can be running using script_epsilon (chmod 777 script_epsilon, then ./script_epsilon)
