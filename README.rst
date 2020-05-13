InvPlane
========

Rotate a planetary system so that an inclination of 0 corresponds to the invariable or fundamental plane.



To compile: 

.. code-block:: bash

  gcc -o invplane inv_plane.c -lm
                                                                   
To execute:

.. code-block:: bash

  invplane [-amhv] <infile>

The infile  contains N+2 lines with the format

.. code-block:: bash

  NumberOrbiters
  CentralMass
  Mass SemimajorAxis Eccentricity Inclination LongitudeAscendingNode ArgumentPericenter MeanAnomaly
  ...
  
where NumberOrbiters is the number of orbiters (don't count the primary!). The next line is the mass of the primary in solar masses. The next N lines contain the orbital elements of the orbiters in the order shown. The units are solar masses, AU, and degrees.                                                       
                                                                   
The user specifies the output format at the command line: 

====   ============
-a     an ASCII text file                                            
-m     a big.in file for use with MERCURY                            
-h     a .hnb file for use with HNBODY                               
-v     report all warnings                                           
====   ============


For example, to transform a system with parameters provided in a file called system.in into a big.in file for use with MERCURY,     
the command is: 

.. code-block:: bash

  invplane -m system.in.                             
                                                                   
Thanks to Russell Deitrick, Pramod Gupta and Joseph Livesey for testing.  
