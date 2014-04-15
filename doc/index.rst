faststat
========

faststat is a *streaming*, *light-weight* statistics library designed for embedding
in other Python applications.  *Steaming* means
that faststat operates on data points as they arrive, without needing to store
previous data points.  *Light-weight* means that statistics do not take up a great
deal of CPU or RAM.  Adding a data point to a stat object is a 0.5 - 3 microsecond operation.  
Each stat object takes about 4kiB of memory.

.. toctree::
   :maxdepth: 2

Basic usage
-----------
The most basic usage of faststat is to create a Stats object to represent some continuous
variable, and then add data points to it.  

::

   >>> import faststat
   >>> s = faststat.Stats()
   >>> for i in range(100):
   ...    s.add(i)
   ...
   >>> s
   <faststat.Stats n=100 mean=49.5 quartiles=(23.4, 49.6, 73.8)>

Stats class
-----------


Examples
--------

Tracking the average number of items present each time a new item was added.

.. code-block:: python

   import collections
   import faststat

   class Queue(object):
      def __init__(self):
         self.deq = collections.deque()
         self.put_stats = faststat.Stats()

      def put(self, item):
         self.put_stats.add(len(self.deq))
         self.deq.append(item)

      def get(self):
         return self.deq.popleft()

