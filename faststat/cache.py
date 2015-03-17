class LRUCache(object):
    '''
    Implements an LRU cache based on a linked list.
    Performance is about 1.1 microseconds per set/get on a core i7
    '''
    def __init__(self, maxlen=10000):
        self.map = {}
        self.root = root = []
        root[:] = [root, root]
        self.maxlen = maxlen
 
    def __getitem__(self, key):
        val = self.map[key][3]
        self[key] = val
        return val
 
    def add(self, key, val):
        # node[0] = prev; node[1] = next
        root = self.root
        if key in self.map:
            # remove from map
            link = self.map.pop(key)
            # remove from list
            link[0][1], link[1][0] = link[1], link[0]
        else:
            link = [None, None, key, val]
        discard = None
        if len(self.map) >= self.maxlen:
            # pop and discard the oldest item
            discard = root[0]
            discard[0][1], root[0] = root, discard[0]
            self.map.pop(discard[2])
        # insert into map
        self.map[key] = link
        # insert into list
        link[0], link[1] = root, root[1]
        root[1][0] = link
        root[1] = link
        if root[0] is root:
            root[0] = root[1]
        if discard:
            return discard[2], discard[3]
 
    __setitem__ = add
 
    _unset = object()
    
    def pop(self, key=_unset):
        # remove from map and list
        if key is LRUCache._unset:
            link = self.root[0]
            self.map.pop(link[2])
        else:
            link = self.map.pop(key)
        link[0][1], link[1][0] = link[1], link[0]
        return link[3]
 
    def __contains__(self, key):
    	return key in self.map
 
    def __len__(self):
        return len(self.map)
 
    def items(self):
        return [(k, self.map[k][3]) for k in self.map]
 
    def keys(self):
        return self.map.keys()
 
    def values(self):
        return [self.map[k][3] for k in self.map]
 
 
class SegmentedCache(object):
    '''
    Implements a Segmented LRU cache based on an LRU cache.
    '''
    def __init__(self, maxlen=10000):
        self.probationary = LRUCache(maxlen)
        self.protected = LRUCache(maxlen / 2)
        self.maxlen = maxlen
 
    def __getitem__(self, key):
        if key in self.protected.map:
            # already protected, nothing to do
            return self.protected[key]
        if key in self.probationary.map:
            # promote to protected
            val = self.probationary.pop(key)
            discard = self.protected.add(key, val)
            if discard:
                self.probationary.add(discard[0], discard[1])
            return val
        raise KeyError(key)
 
    def add(self, key, val):
        if key in self.protected:
            self.protected[key] = val
        elif key in self.probationary:
            self.probationary.pop(key)
            discard = self.protected.add(key, val)
            if discard:
                self.probationary.add(discard[0], discard[1])
        else: # totally brand new key being added
            self.probationary.add(key, val)
            if len(self.probationary.map) + len(self.protected.map) > self.maxlen:
                self.probationary.pop()
 
    __setitem__ = add
 
    def __contains__(self, key):
    	return key in self.probationary or key in self.protected
 
    def __len__(self):
        return len(self.protected) + len(self.probationary)
 
    def items(self):
        return self.protected.items() + self.probationary.items()
 
    def keys(self):
        return self.protected.keys() + self.probationary.keys()
 
    def values(self):
        return self.protected.values() + self.probationary.values()
 
 
if __name__ == "__main__":
	cache_size = 7
	sg = SegmentedCache(cache_size)
	r = range(10000)
	for i in r:
		sg[i] = i
	for i in r[-cache_size:]:
		assert i in sg
 
	import time
 
	r = range(int(5e5))
	s = time.time()
	for i in r:
		sg[i] = i
		sg[i] = i
	print "{0:.2f}us".format(time.time() - s)
