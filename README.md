# QCrop
QCrop for the NYUAD Quantum Computing Hackathon for Social Good in the Arab world

Contains functions for both non error mitiagted and mitigated circuits.

Initialization:
```crop=QCrop()```

Non-error mitigated circuits:
```
crop.generateGraph(molecule,lower_bound,upper_bound,step)
crop.plotGraph()
```

Error Mitigated circuit:
```crop.ErrorMitigatedCircuit(molecule,distance)```
