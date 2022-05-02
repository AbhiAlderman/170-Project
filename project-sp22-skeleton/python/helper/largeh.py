import pickle

d = 100
rp = 14
rs = 3

dictrp = {}
dictrs = {}
for i in range(d): #tower penalty = 8
    for j in range(d):
        lstrp = []
        lstrs = []
        for x in range(d):
            for y in range(d):
                if (i - x)**2 + (y - j)**2 <= rp**2:
                    lstrp.append((x,y))
                if (i - x)**2 + (y - j)**2 <= rs**2:
                    lstrs.append((x,y))
        dictrp[(i, j)] = lstrp
        dictrs[(i, j)] = lstrs
with open('largepenalty.pkl', 'wb') as f:
    pickle.dump(dictrp, f)
with open('largesignal.pkl', 'wb') as f:
    pickle.dump(dictrs, f)