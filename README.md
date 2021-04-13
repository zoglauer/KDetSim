# How to get KDetSim to work in 2021

Clone the repository
```
git clone https://github.com/zoglauer/KDetSim
```

Add the directory to your LD_LIBRARY_PATH. For bash on Linux this would be:
```
export KDETSIM="<path to>/KDetSim"
export LD_LIBRARY_PATH=${KDETSIM}:${LD_LIBRARY_PATH}
```

Run cmake & make
```
cd KDetSim
cmake
make
```

Switch to the examples directory, and run one of the TestStripDetecror examples):
```
root -x TestStripDetector_1.C
```

