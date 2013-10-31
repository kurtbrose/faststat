import faststat

NUM = 40

s = faststat.Stats()

for i in range(NUM):
   s.add(i)

print "0->", NUM, sorted(s.get_percentiles().values())

s = faststat.Stats()

for i in reversed(range(NUM)):
	s.add(i)


print NUM, "->0", sorted(s.get_percentiles().values())
