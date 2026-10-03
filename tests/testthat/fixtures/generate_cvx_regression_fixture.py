"""Independent numerical regression fixture from NeuroKit2, no Physio calls."""
import os, json, hashlib, inspect, platform
from pathlib import Path
p=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(p/'mplconfig')
import numpy as np
import pandas as pd
import scipy
import neurokit2 as nk
import cvxopt
from neurokit2.eda.eda_phasic import _eda_phasic_cvxeda
assert nk.__version__ == '0.2.13'
n=120; sr=4.; tt=np.arange(n)/sr
y=2+.05*np.sin(2*np.pi*tt/20)+.002*tt
for at in [4.,12.,21.]:
    z=np.maximum(tt-at,0)
    y+=.8*(np.exp(-z/2)-np.exp(-z/.7))
y+=.004*np.cos(2*np.pi*tt*1.3)
mu=y.mean(); sd=y.std(ddof=1)
settings=dict(sampling_rate=sr,tau0=2.,tau1=.7,delta_knot=10.,alpha=.0008,gamma=.01,reltol=1e-12)
solver_options=dict(abstol=1e-12,feastol=1e-12)
cvxopt.solvers.options.update(solver_options)
tonic,phasic=_eda_phasic_cvxeda((y-mu)/sd,**settings)
out=p/'cvxeda_neurokit2_fixture.csv'
pd.DataFrame(dict(time_sec=tt,input=y,tonic=tonic*sd+mu,phasic=phasic*sd)).to_csv(out,index=False,float_format='%.17g')
meta=dict(source='NeuroKit2 0.2.13 neurokit2.eda.eda_phasic._eda_phasic_cvxeda',source_url='https://neuropsychology.github.io/NeuroKit/_modules/neurokit2/eda/eda_phasic.html',reference_doi='10.1109/TBME.2015.2474131',python=platform.python_version(),neurokit2=nk.__version__,numpy=np.__version__,scipy=scipy.__version__,cvxopt=cvxopt.__version__,function_sha256=hashlib.sha256(inspect.getsource(_eda_phasic_cvxeda).encode()).hexdigest(),fixture_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),parameters=settings,additional_cvxopt_solver_options=solver_options,samples=n,input_standardization=dict(mean=float(mu),sample_sd=float(sd),ddof=1),output_scaling='tonic = standardised_tonic * sample_sd + mean; phasic = standardised_phasic * sample_sd',physio_parameters=dict(sr=sr,tau1=.7,tau2=2.,delta_knot=10.,alpha=.0008,gamma=.01),max_component_error_tolerance_in_input_sd=1e-4,scope='Numerical implementation agreement only, not physiological recovery accuracy. No decimation. NeuroKit2 retains historical zero-first-two-row ARMA boundary, distinct from author 2025 boundary revision.')
(p/'cvxeda_neurokit2_fixture.json').write_text(json.dumps(meta,indent=2)+'\n')
print(out)
print('input sample SD',sd)
print('dk=2 self-convolved triangle',np.convolve([1,2,1],[1,2,1])/6)
