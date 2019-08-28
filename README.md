# StInLaSimFoam
Hybrid solver between icoFoam and simpleFoam. StInLaSim stands for Steady Incompressible Laminar Simple. 

No need to specify turbulence properties file in the constant folder as we do in SimpleFOAM (OpenFOAMV1906). 

It employs steady state SIMPLE algorithm. 

fvOptions has been enabled which was not present in icoFoam. 
