InvPlane
========

Source code to reference orbital elements to the invariable plane.

To compile: gcc -o invplane inv_plane.c -lm
                                                                   
This code tranforms a system in which inclinations and longitudes of ascending node is measured from an arbitrary reference frame into one in which they are measured from the invariable plane. The input file  contains N+2 lines, where N is the number of orbiters (don't count the primary!). The first line is N, the second line is the mass of the central body, and the next N lines contain the orbital parameters of the orbiters in the format: Mass SemimajorAxis Eccentricity Inclination LongitudeAscendingNode ArgumentPericenter MeanAnomaly. The units are solar masses, AU, and degrees.                                                       
                                                                   
The user specifies the output format at the command line:          
-a   an ASCII text file                                            
-m   a big.in file for use with MERCURY                            
-h   a .hnb file for use with HNBODY                               
-v   report all warnings                                           
                                                                   
For example, to transform a system with parameters provided in a file called system.in into a big.in file for use with MERCURY,     
the command is: invplane -m system.in.                             
                                                                   
Thanks to Russell Deitrick and Pramod Gupta for testing.  
