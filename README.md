faststat
========

fast online statistics collection

very simple API:

```python
>>> import faststat
>>> stats = faststat.Stats()
>>> stats.add(point)
>>> stats.add(point)
...
```


The following properties are accesible on a Stats object: n, min, max, variance, skewness, and kurtosis.
In addition, a Stats object tracks percentiles.

Performance is pretty good: 0.63 microseconds per point on my machine.  (Provided the C module is available.)

In pure Python mode, performance is about 9 microseconds per point.

```python
0.615999937057 microseconds per point
mean (should be 1) 0.998333953189
kurtosis / reference kurtosis -0.0021881144433 -0.00220621681959
variance / reference variance 0.999219190297 0.999219190297
skewness (should be 0) -0.0071960817771
max, min 5.83625092886 -3.4749002526
m2, m3, m4 999218.191078 -7187.64448532 2993126.28574
9.00099992752 microseconds per point
mean (should be 1) 0.998333953189
kurtosis / reference kurtosis -0.0021881144433 -0.00220621681959
variance / reference variance 0.999219190297 0.999219190297
skewness (should be 0) -0.0071960817771
max, min 5.83625092886 -3.4749002526
m2, m3, m4 999218.191078 -7187.64448532 2993126.28574
```
