import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from collections import defaultdict


def pnpe(n, pp, pn, pe):
    if n > 0:
        prob = ((pe**2*pn*pp**2*((n-1)*(pe*(1-pp))**(n-2)))/((pn-pp)*(pe+pp-pe*pp)**(n+2))) + ((pe*pn*pp**2*(pe*(1-pn))**n)/((pn-pp)**2*(pe-pe*pn)*(pe+pn-pe*pn)**(n+1))) - ((pe*pn*pp**2 * (pe*(1-pp))**n * (pe-2*pn + 3*pp + 2* pe*pn - 3 *pe*pp))/((pn - pp)**2*(pe-pe*pp)*(pe + pp - pe*pp)**(n+2)))
        return prob
    else:
        prob = ((pe**2*pn*pp**2*((1/(pe**2*(pp-1)**2))+(((n-1)*((-pe*(pp-1))/(pe+pp-pe*pp))**(n-2))/((pe+pp-pe*pp)**2))))/((pn-pp)*(pe-pp-pe*pp)**2)) - ((pn*pp**2*(pe-1)**3)/((pe+pn - pe*pn)*(pe+pp - pe*pp)**2)) +((pe*pn*pp**2*(((-pe*(pn-1))/(pe+pn-pe*pn))**n-1))/((pn-pp)**2*(pe-pe*pn)*(pe+pn-pe*pn))) - ((pe*pn*pp**2*(((-pe*(pp-1))/(pe+pp-pe*pp))**n-1)*(pe-2*pn+3*pp+2*pe*pn - 3*pe*pp))/((pn - pp)**2*(pe-pe*pp)*(pe+pp-pe*pp)**2))
        return prob

test = pd.read_csv("Nt_seq.csv")
test = test.loc[test['lens_n1'] <=46] #удалены выбросы - длины > 46
x = list(test.iloc[:,7].sort_values())

freq = defaultdict(int)
for l in x:
        freq[l] += 1
freq = dict(freq)
x_exp = list(freq.keys())
y_exp = list(freq.values())

y_exp_norm = [y/sum(y_exp) for y in y_exp] #разделила каждое число встречаний на общее число образцов

y = [pnpe(z,0.2, 0.02715272, 0.15577882) for z in x_exp]

#plot
plt.figure()
plt.plot(x_exp, y_exp_norm, color = "red")
plt.plot(x_exp,y)
plt.grid()
plt.savefig('7_model.png', bbox_inches='tight')

#cord blood

test = pd.read_csv("output_periph.csv")
test = test.loc[test['len'] <=46]
def likelihood(params):
    pp, pn, pe = params[0], params[1], params[2]
    return -np.sum([(int(test.freq.loc[test.len == x0])/sum(test["freq"]))*np.log(pnpe(x0, pp, pn, pe)) for x0 in list(test.len.loc[test.len > 0])])

y = [pnpe(z,0.21308797, 0.08278096, 0.37583493) for z in list(test.len)]

#plot
plt.figure()
plt.plot(list(test.len), list(test.freq/sum(test["freq"])), color = "red")
plt.plot(list(test.len),y)
plt.title("Periphery")
plt.grid()
plt.savefig('bm_plot1.png', bbox_inches='tight')






#likelihood
def likelihood(params):
    pp, pn, pe = params[0], params[1], params[2]
    return -np.sum([(freq[x0]/sum(y_exp))*np.log(pnpe(x0, pp, pn, pe)) for x0 in x_exp])

   
initial_guess = np.array([0.51,0.49,0.5])
result = minimize(likelihood, x0 = initial_guess, bounds = ((0.2, 0.999), (0.01, 0.999), (0.001,0.999)))
print(result.x)