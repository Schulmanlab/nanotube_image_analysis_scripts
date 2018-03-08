#from pymc3 import *
import pymc



x_data = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
y_data = [0.0, 1.5634750934814385, 1.7524718792124605, 2.0091973537076582, 2.6364488552224215, 2.4038551834251547, 3.230498562161169, 2.9100071880960212, 3.27451910383156, 3.1834151125538765]
y_variance = .015

data = dict(x = x_data, y = y_data)

'''with Model() as model: # model specifications in PyMC3 are wrapped in a with-statement
    # Define priors
    #sigma = HalfCauchy('sigma', beta=10, testval=1.)
    sigma = y_variance 
    intercept = 0
    x_coeff = Normal('x', .25, sd=.3)

    # Define likelihood
    likelihood = Normal('y', mu=intercept + x_coeff * x_data,
                        sd=sigma, observed=y_data)

    # Inference!
    trace = sample(3000, cores=2) # draw 3000 posterior samples using
'''



alpha = pymc.Uniform('alpha', lower=0.1, upper=1)

x = pymc.Normal('x', mu=0, tau=1, value=x_data, observed=True)

@pymc.deterministic(plot=False)
def linear_regress(x=x, alpha=alpha):
    return x*alpha

y = pymc.Normal('output', mu=linear_regress, tau = 1.0/y_variance, value=y_data, observed=True)

model = pymc.Model([x, y, alpha])
mcmc = pymc.MCMC(model)
mcmc.sample(iter=100000, burn=10000, thin=10)

print alpha.summary()