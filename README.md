faststat
========

fast online statistics collection

very simple API:

'''python
>>> import faststat
>>> stats = faststat.Stats()
>>> stats.add(point)
>>> stats.add(point)
...
'''


The following properties are accesible on a Stats object: n, min, max, variance, skewness, and kurtosis.

Performance is pretty good: 0.24 microseconds per point on my machine.  (Provided the C module is available.)

In pure Python mode, performance is about 3 microseconds per point.

'''python
>>> import faststat
>>> faststat.test()
0.233999967575 microseconds per point
mean (should be 1) 0.998333953189
kurtosis (should be 0) -2.88310762388
variance (should be 1) 0.999219190297
skewness (should be 0) -0.0071960817771
'''
