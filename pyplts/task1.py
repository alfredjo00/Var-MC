import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

sns.set('talk')

## Load Data

data = np.loadtxt('data/positions.txt')[1500:,:]

E = np.loadtxt('data/energies.txt')

corrs = np.loadtxt('data/corr.txt')[1500:]

def normalize(x, y):
    norm = np.trapz(y, x=x)
    return y/norm

################################################################

def plot_distributions():
    a0 = 5.29177210903e-1 #Angstrom
    x = np.array( [d for d in data[1500:,0] if 0<d] )  * a0
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    #ax.hist(x,bins=100, density=True, label="Sampled $ | \psi_{T} | ^2 $")
    
    rspace = np.linspace(0,4.5,200)
    
    rho = lambda Z: [ Z **3 * 4 * r**2 * np.exp(-2 * Z * r) for r in rspace]
    
    rho_2 = normalize(rspace * a0, rho(2))
    ax.plot(rspace * a0, rho_2, '--k', label=r'$ \rho(r),  Z=2$')
    
    rho_27 = normalize(rspace * a0, rho(27/16))
    ax.plot(rspace * a0, rho_27, '-.r',  label=r'$ \rho(r), Z=27/16$')
    
    
    y, bins= np.histogram(x, bins=500, density=True)    
        
    y_p = [p * r **2 for p,r in zip(y,bins[1:])]
    y_p_norm = normalize(bins[1:], y_p)
    ax.plot(bins[1:],y_p_norm, label="$| \psi_{T}(r) | ^2 r^2 $ (normalized)")
    
    ax.legend()
    
    
    ax.set_xlabel(u'Distance $r$ from origo (\u00C5)')
    ax.set_ylabel(r"Probability density")
    
    fig.tight_layout()
    fig.savefig("./graphs/task1_test_dist.jpg", dpi=150)

################################################################

def correlations():

    fig,ax = plt.subplots(figsize=(8,6))
    ax.hist(corrs, bins=35, density=True, label=r'Sampled $|\Psi_{T}(\mathbf{r}_1, \mathbf{r}_2)|^2$')
    ax.plot([-1,1], [0.5,0.5], '--k', label='Uncorrelated data')
    ax.set_ylabel('Probability density')
    ax.set_xlabel(r"$\cos\theta$")
    ax.legend()
    fig.tight_layout()
    
    fig.savefig('graphs/corrs.jpg')

################################################################

def local_energies():
    def moving_average(a, n=150):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    
    fig,ax = plt.subplots(figsize=(8,6))
    ax.plot(E[:3000], label='$E_L$')
    ax.plot(moving_average(E[:3000]), '-k', linewidth=5, label="Moving average of $E_L$")
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r"Local energy $E_L$ (Hartree $E_\mathrm{h}$)")
    ax.legend()
    fig.tight_layout()
    fig.savefig('graphs/E_eq.jpg')

################################################################

def alpha_sweep():
    alphas = np.loadtxt("data/alpha_sweep.txt")
    
    fig, ax  = plt.subplots(figsize=(8,6))
    x, y, yerr = alphas[:,0], alphas[:,1], alphas[:,2]
    ax.set_title(r'Local energy average, varying $\alpha$')
    ax.errorbar(x, y, 2 * yerr, fmt='--b', ecolor='k', capsize=7,
                 linewidth=2, label=r'$\mu \pm 2\cdot \sigma$',
                 errorevery=3)
    
    ax.set_xlabel("$\alpha$")
    ax.set_ylabel("$\langle E_L  \rangle$")
    
    x_fit = x[2:30]
    y_fit = y[2:30]
    
    z = np.polyfit(x_fit, y_fit, 2)
    p = np.poly1d(z)
    
    ax.plot(x_fit, p(x_fit), '-.r', label=f'Quadratic fit, min = {x[np.argmin(p(x_fit))]:.5f}', linewidth=2.5)
    #errorbar(x, y, yerr, xerr, fmt, ecolor, elinewidth, capsize, barsabove, lolims, uplims, xlolims, xuplims, errorevery, capthick)
    
    
    
    ax.set_ylabel(r"$E_L(\mathbf{r}_1, \mathbf{r}_2)$ (Hartree $E_\mathrm{h}$)")
    ax.set_xlabel(r"Parameter $\alpha$")
    ax.legend()
    fig.tight_layout()
    fig.savefig("graphs/alpha_sweep.jpg")

################################################################

def grad_descent_plot():
    
    grad_data = np.loadtxt('data/grad_descent.txt')
    beta0_01 = np.loadtxt('data/grad_descent_0_01.txt')
    beta0_1 = np.loadtxt('data/grad_descent_0_10.txt')
    
    x = grad_data
    
    alpha_mn = np.mean(x[-15:])
    alpha_std = np.std(x[-15:])
    print(alpha_mn,alpha_std)
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(x, '-.ok', markerfacecolor='none', label=r'Gradient descent, $\beta = 0.75$')
    
    ax.plot([0, len(x)-1], [alpha_mn, alpha_mn], '--b', label=fr'Optimal $\alpha = {alpha_mn:.5f} \pm {alpha_std:.5f} $')
    
    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'$\alpha$')
    ax.legend()
    ax.set_ylim(0.14,0.162)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))    
    fig.tight_layout()
    fig.savefig('graphs/grad_descent.jpg')
    
################################################################

def plot_corr_fn():
    
    fn_data = np.loadtxt('data/correlation_fn.txt')
    
    y = fn_data    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(np.arange(len(y)), fn_data, '-k', markerfacecolor='none', label=r'$\Phi(k)$')
    
    ind_arg = np.argmin(abs( y - np.exp(-2)))       
    
    ax.plot(np.arange(len(y)), np.ones(len(y)) * np.exp(-2), '--b', label=fr'$\Phi(k) = e^{{-2}}, k \approx {ind_arg}$' )
    ax.plot([ind_arg, ind_arg], [0, np.exp(-2)], '--b')
    
    
    ax.set_title(r'Correlation function')
    ax.legend()
    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'Correlation')
    
    fig.tight_layout()
    fig.savefig('graphs/correlation_fn.jpg')
    
################################################################


def plot_blocks():
    
    block_data = np.loadtxt('data/blocks.txt')
    
    y = block_data
    
    block_mn = np.mean(y[-20:])
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(np.arange(len(y)), y, '-k', markerfacecolor='none', label='$n_s$')    
    
    ax.plot(np.arange(len(y)), np.ones(len(y)) * block_mn, '--b', label=fr'Converged $n_s = {block_mn:.3f}$')
    
    ax.legend()
    
    ax.set_title(r'Block averaging')
    
    ax.set_xlabel('$M_B$')
    ax.set_ylabel(r'Statistical inefficiency')
    
    fig.tight_layout()
    fig.savefig('graphs/blocks.jpg')
    
################################################################

if __name__ == '__main__':
    #alpha_sweep()
    local_energies()
    correlations()
    plot_distributions()
    grad_descent_plot()
    plot_corr_fn()
    plot_blocks()


