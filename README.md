# Nonparametric and Neural Nets regression
Some neural networks implementations in Python and R from my Master's Thesis.

To test the Python algo (fastest) you need python3 (recommended in a virtualenv) and once you have done that just do
```
source YourVirtualEnvFolder/bin/activate
pip install -r requirements.txt
```

In order to run an example you do
```python test.py```

More examples are commented out below the one which will run by default.

The R implementation is a general R project which you can build and write a custom script to test it. It contains nonparametric kernel regression algorithms as well. The neural networks implementations is for just one and two hidden layers and not as general as the python implementation. If you build the project the roxygen files which guide you through what to do will be generated as well.

The nonparametric iv regression is the sieve implementation of Horowitz and Lee (2007) and Horowitz(2012) with the uniform confidence bands under both the monotonicity and Lipschitz assumption on the unknown g. It does as well the MCMC simulation for the testing of the theory. You just fire it up and it runs the MCMC.
