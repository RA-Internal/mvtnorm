Revolution Analytics package "mtvnorm" validation. 

--Runit tests constructed in scripts "testRunScript.R" and "RU.pmvnorm.R"
found in the "tests" folder. (run script, "testRunScript.R" with folder as workspace)
--To compare runs under versions 3.0.1 and 3.0.2, check values listed in files:
mtvnorm301.Rout
mtvnorm302.Rout
	
--Two "failed" tests:  When generating normally distributed random numbers
and setting the seed for future exact replication, be sure to use either the Cholesky 
decomposition, or the other two decompositions.  cor(x1,x2) < 1 (roughly .985), if
one vector x is produced using a Cholesky decomposition and the next vector is not.

# note:  documentation for mvtnorm differs mildly at points from CRAN guidelines
# note:  The large list of reverse depends strongly suggests this package functions well.
