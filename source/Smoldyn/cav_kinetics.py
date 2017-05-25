import numpy as np
v=np.linspace(-100,100,500)

inf_n=1.0/(1.0+np.exp((v+15.0)/-8.0))
tau_n=0.06+0.75*4*np.sqrt(0.55*0.45)/(np.exp((v+15.0)*0.55/8.0) + np.exp(-(v+15.0)*0.45/8.0));
rate_open=inf_n/tau_n
rate_close=(1-inf_n)/tau_n
