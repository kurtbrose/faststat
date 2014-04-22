import collections

# This is a fairly naive implementation of the algorithm described in
# "Effective Computation of Biased Quantiles over Data Streams"

class Quantiles(object):
    def __init__(self, f=None):
        self.points = None
        self.n = 0
        self.f = f or biased_quantiles_f

    def insert(self, val):
        if self.points is None or val < self.points.val:
            # less than first
            self.points = _point(val, 1, 0, self.points)
            return
        prev = self.points
        cur = self.points.next
        r = 0
        while cur and cur.next:
            if val < cur.val:
                break
            r += cur.delta
            prev = cur
            cur = cur.next
        else:  # ran off the end
            (cur or prev).next = _point(val, 1, 0, None)
            return
        new = _point(val, 1, max(int(self.f(r, self.n)) - 1, 0), cur)
        prev.next = new
        self.n += 1
        if self.n % COMPRESS_INTERVAL == 0:
            self.compress()

    def compress(self):
        pointlist = []
        cur = self.points
        r = [0]
        while cur:
            pointlist.append(cur)
            r.append(r[-1] + cur.delta)
            cur = cur.next
        for i in range(len(pointlist) - 2, 0, -1):
            if pointlist[i].delta + pointlist[i + 1].delta + pointlist[i + 1].width <= self.f(r[i], self.n):
                # merge
                pointlist[i].next = pointlist[i + 1].next
                pointlist[i].val = pointlist[i + 1].val
                pointlist[i].delta += pointlist[i + 1].delta
                pointlist[i].width = pointlist[i + 1].width

    def query(self, p):
        cur_r = 0
        cur_point = self.points
        while cur_point:
            cur_r += cur_point.delta
            if (cur_r + cur_point.next.delta + cur_point.next.width > 
                        self.n * p + self.f(self.n * p, self.n) / 2):
                return cur_point, (cur_r + cur_point.delta, cur_r + cur_point.delta + cur_point.width)
            cur_point = cur_point.next

    def get_pointlist(self):
        if not self.points:
            return []
        p = [self.points]
        while p[-1].next:
            p.append(p[-1].next)
        return p



class _point(object):
    def __init__(self, val, delta, width, next):
        self.val, self.delta, self.width, self.next = val, delta, width, next

    def __repr__(self):
        return "_point(val={0}, delta={1}, width={2})".format(self.val, self.delta, self.width)


def biased_quantiles_f(r_i, n):
    return 2 * ERROR_RATE * r_i


def targeted_quantiles_f(percentiles_errors):
    def f(r_i, n):
        bounds = []
        for p, e in percentiles_errors:
            if r_i < p * n:
                bounds.append(2 * e * (n - r_i) / (1 - e))
            else:
                bounds.append(2 * e * r_i / p)
        return min(bounds)
    return f


ERROR_RATE = 0.001
COMPRESS_INTERVAL = int(1 / ERROR_RATE)

# val is the current (approximate) value

# delta is the difference between the lowest possible rank of the current
# point/value and the previous point

# width is the differencet between the lowest and highest possible rank
# of the current point/value

# this data structure ensures that new points can be inserted into
# the middle of the linked list


# performance of naive algorithm is very bad -- 300 - 700 microseconds
# (0.3 to 0.7 ms).  this is about 20-40x slower than python piece-wise
# parabolic algorithm;  ~300x slower than C piece-wise parabolic
def test(q=None):
    import random, time

    data = [random.normalvariate(1.0, 1.0) for i in range(int(1e4))]
    q = q or Quantiles()
    start = time.time()
    for d in data:
        q.insert(d)
    print (time.time() - start) * 1e6 / len(data), "microseconds per point"
    return q

# about 400 microseconds per point
def test_targeted():
    TARGETS = ((0.25, 0.001), (0.5, 0.001), (0.75, 0.001), (0.9, 0.001), (0.95, 0.001), (0.99, 0.001))
    f = targeted_quantiles_f(TARGETS)
    return test(q=Quantiles(f))

