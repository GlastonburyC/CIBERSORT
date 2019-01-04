# Using CIBERSORT to estimate cell type composition.


### CIBERSORT used to estimate fractions:

``` 
R –-no-restore
library(Rserve)
Rserve()
q()
java -jar CIBERSORT.jar -M TwinsUK.adiposeSamples.txt -B Adipose.sigMatrix.txt –n 1000 >> TwinsUK.Adipose.CellEsts.txt
```
