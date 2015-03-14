# Single Input Files #

shab-pimd accepts an input file `InputFileName.in`.

example.in:
```
M = 32
2 1 32 10.0 0.01 10.0 100 1000 10000 0 2 0 3 4 5
2 1 32 10.0 0.01 10.0 100 1000 10000 1 2 0 3 4 5
2 1 32 10.0 0.01 10.0 100 1000 10000 2 2 0 3 4 5
```

Here the first line is the input file label:
```
M = 32
```
implying there are `32` Beads as can be seen in the next few lines which contain the settings for three different simulations. This label will be used later when making plots.

These settings are:
```
nPart nD nBead beta dt L nBins eSteps rSteps transformation thermostat interaction nNH SYOrder nNHsteps
```

We can pass these settings to shab-pimd in the following way:
```
./shab-pimd inputs/example.in 3
```
where the first argument is the input file location, and the second argument in which settings line to use. Note that this means `3` actually uses the fourth line of `example.in` since the first line is the label.

Each run puts all traces in `data/traces/`.

# Batch Runs #

If you wish to run every line from an input file, you may use the `BatchRun.sh` BASH script as follows:
```
sh BatchRun.sh 2 example.in
```
where `2` denotes the number of processes to run simultaneously (if your computer has more than one processor), and `example.in` is the input file name. Note that `BatchRun.sh` automatically checks the `inputs/` folder for the input file.

`BatchRun.sh` can also handle multiple input files:
```
sh BatchRun.sh 2 example.in example.in example.in
```

The to screen output from each simulation is found in `data/output/`, while again each run puts all traces in `data/traces/`.

# Plotting Data #

The BASH script `AllPlots.sh` also accepts input files as follows:
```
sh AllPlots.sh Example example.in
```
where `Example` is the label of the plot file, and `example.in` is the name of the input file in `inputs/`. This script generates a single plot for each observable measured. These figures are stored in `data/figures/`.

**WARNING:** Right now this assumes the x-axis is beta ranging from 1 to the number of data points. You will want to change the section under _# Beta Values_ if this is not the case.

`AllPlots.sh` can also accept multiple input files:
```
sh AllPlots.sh Example example.in example.in
```
Again, only one figure per observable will be created, thus in this example each plot will have two lines. The legend entry for each line is the input file label mentioned above.