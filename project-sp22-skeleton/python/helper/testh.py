import pickle

d = 10
rp = 2
rs = 1

dictrp = {}
dictrs = {}
for i in range(d): #tower penalty = 8
    for j in range(d):
        lstrp = []
        lstrs = []
        for x in range(d):
            for y in range(d):
                if (i - x)**2 + (y - j)**2 <= rp**2 and not (i==x and j==y):
                    lstrp.append((x,y))
                if (i - x)**2 + (y - j)**2 <= rs**2:
                    lstrs.append((x,y))
        dictrp[(i, j)] = lstrp
        dictrs[(i, j)] = lstrs
with open('testpenalty.pkl', 'wb') as f:
    pickle.dump(dictrp, f)
with open('testsignal.pkl', 'wb') as f:
    pickle.dump(dictrs, f)