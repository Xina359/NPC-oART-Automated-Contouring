# Terminology

Use **automated contouring** for the clinical system, its intended radiotherapy
task, the repository title, workflow descriptions and the algorithm-access
statement. Use **contours** for structures reviewed and edited by clinicians.
For target-volume definition, **GTV delineation** and **CTV delineation** are
appropriate task descriptions.

Use **segmentation** for voxel-level prediction, segmentation networks, masks,
training objectives and evaluation functions. Technical identifiers such as
`segmentation_metrics.py` remain unchanged. Use **geometric agreement** when
describing what DSC, HD95 and ASD quantify between automatic and reference
contours.

These are contextual choices, not a claim that contouring is two-dimensional
or that segmentation is exclusively three-dimensional. Both terms are used
in radiotherapy; avoid switching terms merely for stylistic variation within
one description of the same task.
