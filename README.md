# SimplexAlgorithm
 A simplex optimizer written in Python. Final Project for CS 257.
 
Download test files (using Tiny-C compiler), starting from project directory. 

```console
$ cd mps_decomp
$ curl -O https://netlib.org/lp/data/emps.c
$ tcc -c emps.c
$ tcc -o emps emps.o
# cd ../data/netlib
$ curl -O "https://netlib.org/lp/data/{25fv47,adlittle,afiro,agg,agg2,agg3,bandm,beaconfd,blend,boeing1,boeing2,bore3d,brandy,capri,e226,forplan,grow15,grow22,grow7,israel,kb2,lotfi}"
$ for f in 25fv47 adlittle afiro agg agg2 agg3 bandm beaconfd blend boeing1 boeing2 bore3d brandy capri e226 forplan grow15 grow22 grow7 israel kb2 lotfi; do ./../../mps_decomp/emps "${f}" > "${f}.mps"; 
done
```
