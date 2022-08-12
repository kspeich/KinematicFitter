# KinematicFitter
Kinematic Fitting algorithms for a CMS $H\rightarrow aa\rightarrow b\overline{b}\tau^-\tau^+$ analysis

By Kodai Speich, Princeton University

This mostly-C++ repository builds off the kinematic fitting algorithm that already exists in ```/cmssw/PhysicsTools/KinFitter```, and can be cloned into ```/cmssw/PhysicsTools```.

This repository contains three directories.  ```interface``` contains the header files for all of the classes, ```src``` contains the corresponding source files, and ```test``` contains C++ macros as well as histograms and plots.  The classes contained in this repository are mostly particle and constraint classes -- derived classes of ```TAbsFitParticle``` and ```TAbsFitConstraint```, respectively.  These particles and constraints can be added to a ```TKinFitter``` object and then fit.

The ```KinFitOutputModule``` class reads a TTree branches by event, and creates a list of events, which contains the particle-level data.  It then runs the kinematic fitting algorithm on each event and populates a fitted events list.  It can then use these two lists to create histograms.  ```SVFitKinFitOutputModule``` is derived from ```KinFitOutputModule```, and imposes different constraints onto the SVFit values.

The ```KinFitEfficiency``` class reads a ```KinFitOutputModule``` object and calculates $S/B$ and $S/\sqrt{S+B}$ ratios for the $\tau^-\tau^+$ and $b\overline{b}$ reconstruction of the pseudoscalar mass $m_a=45$ GeV as well as the all-particle reconstruction of the Higgs mass $m_h=125$ GeV.

Each object is called in runKinFitter.cc (located in ```test```), which can be run with ```root runKinFitter.cc``` in ```test```.  Another macro, ```plotOutput.cc``` overlays plots of the histograms created in ```KinFitOutputModule```, and can be run by ```root plotOutput.cc``` in ```test```.  To compile the code, run ```scram b -j 8``` in ```KinematicFitter```.
