# SimplexAlgorithm
 A simplex optimizer written in Python. Final Project for CS 257.
 
Download test files (using Tiny-C compiler): 

```console
$ curl -O https://netlib.org/lp/data/emps.c
$ tcc -c emps.c
$ tcc -o emps emps.o
$ curl -O "https://netlib.org/lp/data/{25fv47,80bau3b,adlittle,afiro,agg,agg2,agg3,bandm,beaconfd,blend,bnl1,bnl2,boeing1,boeing2,bore3d,brandy,capri,cycle,czprob,d2q06c,d6cube,degen2,degen3,dfl001,e226,etamacro,fffff800,finnis,fit1d,fit1p,fit2d,fit2p,forplan,ganges,gfrd-pnc,greenbea,greenbeb,grow15,grow22,grow7}"
$ for f in 25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2 bore3d brandy capri cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7; do ./../../mps_decomp/emps "${f}" > "${f}.mps"; 
done
```
