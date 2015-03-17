import math

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


import math


class DistributedQuantiles(object):
    def __init__(self, error):
        self.size = math.ceil(1 / error)
        self.inbuf = []
        self.qs = []
        self.n = 0

    def insert(self, val):
        self.inbuf.append(val)
        if len(self.inbuf) >= self.size:
            self.inbuf.sort()
            i, j = 0, 0
            new_qs = []
            if self.inbuf[0] < self.qs[0][0]:
                while self.inbuf[i] < self.qs[j][0]:
                    new_qs.append((self.inbuf[i], i, i))
            while i < len(self.inbuf) and j < len(self.qs):
                nxtmin = i + self.qs[j][1]
                nxtmax = i + self.qs[j][2] - 1
                if self.inbuf[i] < self.qs[j][0]:
                    nxtval = self.inbuf[i]
                    i += 1
                elif self.qs[j] < self.inbuf[i]:
                    nxtval = self.qs[j]
                    j += 1
                else:
                    val = self.inbuf[i]
                    if len(inbuf) > i + 1 and inbuf[i + 1] == val:
                        if len(self.qs) > j + 1 and self.qs[j + 1] == val:
                            new_qs.append((val, ))
                    # the paper does not cover this case so a supplemental
                    # proof is provided
                    pass
                new_qs.append((nxtval, nxtmin, nxtmax))
            self.qs = new_qs
            self.prune()

    def prune(self):
        pass


'''
Supplemental proof of correctness of equal value element case.
(notation is kept as close as possible to the original paper
within the limitations of ASCII)

MOTIVATION:
The error bound as given in the paper is dependent that
when merging any Q' and Q'', the next-lowest y_s of Q'' and
next-highest y_t of Q'' to any x of Q' must be consecutive.

However, in the general case there may be many elements whose
value are equal.  Consider
Q' = [...,-1,0,1,...]
Q'' = [...,-1,0,0,1,...]

For x of Q' = 0, y_s = -1 and y_t = 1.
y_s and y_t are not consecutive.
The merge operation is no longer correct, because inductively
rmax_Q''(y_s) - rmin_Q''(y_t) <= 3 eps n'' ; which is not <= 2 eps n''

We therefore define a new merge operation for elements of Q' and Q''
which are equal.

This merge operation will also inductively guarantee that only
two consecutive elements may share the same value.


INTUITION:
The intuitive explanation of the meaning of multiple elements
with the same value is that it represents an unbroken
sequence of observations with that same value.  The first
element represents the rank of the beginning of this sequence.
The second element represents the rank of the end of this
sequence.

Because of this, there should never be more than two elements
with the same value.  (If there is a third element, it can be
discarded keeping only the maximum and minimum.)

This also means that if there are two elements with the same
value, rmin may be set to rmax for the lower ranked element,
and rmax may be set to rmin for the higher ranked element.

Because every rank between the two elements of the same
value is known to also be of that value, the top-most rank
in the range of the first value and bottom-most rank in the
range of the second value are known to be the exact locations
of an element of the given value.


NEW MERGE OPERATION:

The proposed new merge operation considers 4 points at a time:

x_r, x_r+1, y_t, y_t+1

where x_r is the minimum un-merged element of Q'
      x_r+1 is the next consecutive element of Q'
      y_t is the minimum un-merged element of Q''
      y_t+1 is the next consecutive element of Q''

in order to generate z_i of Q

all points whose value are equal are consumed in a single step

there are four cases:

for notation convenience,
define the ADD operation of x of Q' and y of Q''
be z of Q:
    rmax_Q(z) = rmax_Q'(x) + rmax_Q''(y)
    rmin_Q(z) = rmin_Q'(x) + rmin_Q''(y)

CASE I:
    x_r != y_t:
        apply the merge operation as in the original paper
CASE II:
    x_r = y_t, and x_r != x_r+1, and y_t != y_t+1:
        z_i = ADD(x_r, y_t)
        move on to x_r+1, y_t+1
case III:
    x_r = y_t = x_r+1 != y_t+1:
        generate two points, z_i and z_i'
        z_i = ADD(x_r, y_t)
        z_i' = ADD(x_r+1, y_t)
        apply SHRINK to the pair of points,
        append the result (either one or two points) to Q,
        move on to x_r+2 and y_t+1
case IV:
    x_r = y_t = y_t+1 != x_r+1:
        generate two points z_i and z_i'
        z_i = ADD(x_r, y_t)
        z_i' = ADD(x_r, y_t+1)
        apply SHRINK to the pair of points,
        append the result (either one or two points) to Q,
        move on to x_r+1, y_t+2
case V:
    x_r = y_t = x_r+1 = y_t+1:
        this means that there are two runs of the same value
        in each of the merging sequences.

        generate two points z_i, z_i':
        z_i = ADD(x_r, y_t)
        z_i' = ADD(x_r+1, y_t+1)
        append both to Q
        move on to x_r+2, y_t+2


PROOF:
Let eps be the error bound.
Let Q' and Q'' be two quantile summaries being merged, each of
which inductively have error <= eps.
Let x be an element of Q' and y be an element of Q'' such that.
Let n' be the number of observations covered by Q' and
n'' be the number of observations covered by Q''.

The proposed merge operation is:
rmax_Q(z_i) = rmax_Q'(x) + rmax_Q''(y)
rmin_Q(z_i) = rmin_Q'(x) + rmin_Q''(y)

We wish to show that rmax_Q(z_i+1) - rmin_Q(z_i) <= 2 eps (n' + n'')

let r be the index of x in Q'
let s by the index of y in Q''

Several things must be proven about the new MERGE:
1- that z_i and z_i+1 are within the error bounds when z_i is the 
        result of CASE II, II, IV, or V, and z_i+1 is the result of CASE I
2- that z_i and z_i+1 are within the error bounds when z_i is the
        result of CASE I, and z_i+1 is the result of CASE II, III, IV, or V
3- that z_i and z_i+1 are within the error bound

There are 3 cases:

case I: 

case I: z_i+1 was merged using the new merge operation
case II: z_i+1 was merged using the merge operation defined in the paper
case IIA: z_i+1 came from x_r+1 of Q'
case IIAi: z_i+1 = z_i
case IIAii: z_i+1 > z_i
case IIB: z_i+1 came from y_s+1 of Q''
case IIBi: z_i+1 = z_i
case IIBii: z_i+1 > z_i

(Note, if x_r+1 = x_r and y_s+1 = y_s, this is back to case I)

case I:
    rmax_Q(z_i+1) - rmin_Q(z_i) = rmax_Q'(x_r+1) + rmax_Q''(y_s+1)
                                   - rmin_Q'(x_r) - rmin_Q''(y_s)
                                = rmax_Q'(x_r+1) - rmin_Q'(x_r)
                                   + rmax_Q''(y_s+1) - rmin_Q''(y_s)
    by inductive property of Q' and Q''
                                <= 2 n' eps + 2 n'' eps
                                <= 2 eps (n' + n'')

case IIAi:
    if x_r = x_r+1:
        by consecutive element contraction lemma:
            if x_r and x_r+1 overlap, combine them
            and there is no z_i, only z_i+1

            otherwise, rmin_S(x_r) = rmax_S(x_r)
            and rmin_S(x_r+1) = rmax_S(x_r+1)

            rmax_Q(z_i+1) - rmin_Q(z_i)
            = 

        

   


consecutive elment contraction lemma:
    let s and s' be consecutive elements of summary S,
    such that s = s'
    CASE I:
    if rmin_S(s') < rmax_S(s):
        maximum rank of s' = rmax_S(s)
        minimum rank of s = rmin_S(s')
        [s and s' are now identical, discard one]
    CASE II:
    if rmin_S(s') >= rmax_S(s):
        maximum rank of s' = rmin_S(s')
        minimum rank of s = rmax_S(s)

    let '[' represent rmin_S, and ']' represent rmax_S:
        [  x_r  ]
             [  x_r+1 ]
    let a = the value of x_r and x_r+1;
    We know that somewhere between rmax_S(x_r) and rmin_S(x_r+1)
    there is a sequence of consecutive 'a' values in the data
    set.  Therefore, there is at least one a in that interval.



-----------------------
Notation:
for readability, everything after a "_" should be read as a subscript
For example, z_i+1 should be read z sub(i+1), not z sub(i) + 1

x of Q, should be read "x element of Q"

Also, the original paper assumes that every element is unique.
Therefore, rmin(v) is a function that takes value v and returns
the minimum rank.  Because we allow two elements to have the same
value, v must be considered to be a tuple.  rmin(v) is selecting
one item from that tuple, rmax(v) is selecting another element.

This does not affect the correctness of any of the proofs from the
original paper.

'''



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
COMPRESS_INTERVAL = 10 # int(1 / ERROR_RATE)

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

if __name__ == "__main__":
    print 'biased quantile condition'
    test()
    print 'targeted quantile condition'
    test_targeted()
