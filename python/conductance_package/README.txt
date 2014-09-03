1) To submit EVB simulations with energy sink

Go to directory Prod/

Directory Prod/sub contains script generate.py. Running the script
asks a series of questions to make sure required files are in the
appropriate directories. After running the script:
 - directories Prod/L-?.??/curr/ populates with required input files
   to run the simulations
 - directory Prod/sub populates with scripts prod.sh and 0001_0100.sub
   to 1901_2000.sub. To send the jobs to the queue, run ./prod.sh


2) To run the conductance analysis

Files needed:

(1) topology file
(2) coordinate trajectory
(3) velocity trajectory
(4) PDB file
(5) file containing one-number per atom specifying the cold-atoms

NOTE: files (1), (2), (3) must be of the desolvated and de-ionized
system (no water, no ions)

To test the package, run:
./conductance.py -a test/dhp.top -b test/dhp.00001_03100.crd -c test/dhp.00001_03100.vel -d test/dhp.uns.premin.pdb -e test/indexes.dat -f ./tmp
