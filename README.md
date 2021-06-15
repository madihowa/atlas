# ClusterCalibration
Cluster Calibration for the Atlas Detector

This project was meant to recalibrate the deposited energy in the atlas calorimeters.  

The orignal data used for this project is currently stored [Here](https://www.dropbox.com/s/93xk0vl6ou193lm/Pi.tar?dl=0).

To run this data

```bash
git clone ssh://git@gitlab.cern.ch:7999/nyoung/ml-calorimeter-reweighting.git
cd ClusterCalibration
mkdir EMDataSets
cd EMDataSets
wget https://www.dropbox.com/s/fp8gb2spyce1x1v/Pi.tar?dl=0
tar -xvf Pi.tar?dl=0.1
rm Pi.tar?dl=0.1
cd ../
./run.sh
```

When the code is finished running it will put all the results in a new folder.

If you have any questions don't hesitate to contact me at "nyoung@uoregon.edu"
