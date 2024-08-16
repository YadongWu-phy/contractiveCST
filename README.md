# ContractiveCST
Numerical results on the contractive unitary for Classical Shadow Tomography.
The codes are written in `Julia` programming language and make use of the `QuantumClifford.jl` package.

## Calculation Scripts
[ContractSample.jl](ContractSample.jl) and [CliffordSample.jl](CliffordSample.jl) are utilized to simulate classical shadow tomography protocols with different unitary ensembles, with knowledge of the location where Pauli operators are applied.

[ContractSlid.jl](ContractSlid.jl) and [CliffordSlid.jl](CliffordSlid.jl) are utilized to simulate classical shadow tomography protocols with different unitary ensembles, without knowledge of the location where Pauli operators are applied.

## Data
Data required to reproduce figures in "Contractive Unitary and Classical Shadow Tomography":

`SizeDistrib.mat` illustrates the difference in size distribution between random Clifford ensemble and contractive unitary ensemble.`sample1.mat` displays the numerical result of the classical shadow protocol for size-k Pauli operators.`sampleSlide.mat` exhibits the numerical result of the classical shadow protocol using the 'sliding trick'.
