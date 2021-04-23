%run EBMA and get motion vectors
[u1,v1] = ExhaustiveSearch(8,13);

%run TSS and get motion vectors
[u2,v2] = ThreeStepSearch(8,13);

%run NTSS and get motion vectors
[u3,v3] = NewThreeStepSearch(8,13);
