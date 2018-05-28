import os

fp = open('script.sh', 'w')
d = '/home/jingweih/Desktop/cars'

t1 = 0
t2 = 0
t3 = 0
t4 = 0
t11 = 0
t12 = 0
t13 = 0
t14 = 0
for o in os.listdir(d):
	folder = os.path.join(d,o)
	if os.path.isdir(folder):
		if (o == 'IM-right'):			
			for f in os.listdir(folder):
				if os.path.isfile(os.path.join(folder, f)):
					if f.endswith('.txt'):
						f1 = open(folder + '/'+ f, 'r')
						a1,b1 = map(int, f1.readline().split())
						c1 = map(float,f1.readline().split())[0]
						d1 = map(float,f1.readline().split())[0]
						f2 = open('/home/jingweih/Desktop/cars/QF/' + f, 'r')
						a2,b2 = map(int, f2.readline().split())
						c2 = map(float,f2.readline().split())[0]
						d2 = map(float,f2.readline().split())[0]
						if a1 < a2:
							t1 += 1
						else:
							t11 += 1
						if b1 < b2:
							t2 += 1
						else:
							t12 += 1
						if c1 < c2:
							t3 += 1
						else:
							t13 += 1
						if d1 < d2:
							t4 += 1
						else:
							t14 += 1
fp.close()

print t1,t2,t3,t4
print t11,t12,t13,t14
