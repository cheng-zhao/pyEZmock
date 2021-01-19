#!/usr/bin/env python3
# pyEZmock.py: this file is part of the pyEZmock package.
#
# pyEZmock: Python wrapper of Effective Zel'dovich approximation mock (EZmock).
#
# Github repository:
#       https://github.com/cheng-zhao/pyEZmock
#
# Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import os
from warnings import warn

class pyEZmock:

  def __init__(self, workdir=None,
      exe='/global/u2/z/zhaoc/work/pyEZmock/bin/EZmock.sh',
      pk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/POWSPEC.sh',
      xi_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/FCFC_2PT_BOX.sh',
      bk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/BISPEC_BOX.sh'):
    """
    Initialise the `pyEZmock` class.

    Parameters
    ----------
    exe: str, optional
        Location of the EZmock executable.
    pk_exe: str, optional
        Location of the powspec executable.
    xi_exe: str, optional
        Location of the FCFC executable for periodic boxes.
    bk_exe: str, optional
        Location of the bispec executable.
    """
    if workdir is None:
      raise ValueError('Working directory should be set via `workdir`.')
    if not os.path.isdir(workdir):
      try: os.mkdir(workdir)
      except: raise IOError(f'Cannot create working directory: {workdir}')
    self.workdir = os.path.abspath(workdir)

    if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
      raise IOError(f'Invalid EZmock executable: {exe}')
    self.ez_exe = os.path.abspath(exe)

    if not (os.path.isfile(pk_exe) and os.access(pk_exe, os.X_OK)):
      warn(f'Invalid power spectrum executable: {exe}')
      self.pk_exe = None
    else: self.pk_exe = os.path.abspath(pk_exe)

    if not (os.path.isfile(xi_exe) and os.access(xi_exe, os.X_OK)):
      warn(f'Invalid 2PCF executable: {exe}')
      self.xi_exe = None
    else: self.xi_exe = os.path.abspath(xi_exe)

    if not (os.path.isfile(bk_exe) and os.access(bk_exe, os.X_OK)):
      warn(f'Invalid bispectrum executable: {exe}')
      self.bk_exe = None
    else: self.bk_exe = os.path.abspath(bk_exe)

    # Initialise EZmock parameters
    self.__param = dict(
      boxsize = None,
      num_grid = None,
      redshift = None,
      num_tracer = None,
      pdf_base = None,
      dens_scat = None,
      rand_motion = None,
      dens_cut = None,
      seed = 1,
      omega_m = 0.307115,
      init_pk = '/global/u2/z/zhaoc/work/pyEZmock/data/PlanckDM.linear.pk'
    )

    # Initialise clustering settings
    self.__pk = False           # indicate whether to compute power spectra
    self.__pkl = []             # power spectrum multipoles to be evaluated
    self.__pkgrid = 256         # grid size for power spectra evaluation
    self.__kmax = 0.3           # maximum k for the power spectra
    self.__dk = 0.01            # bin size of k for the power spectra
    self.__pk_ref = None        # reference power spectrum
    self.__pk_col = []          # columns of the reference P(k) multipoles

    self.__xi = False           # indicate whether to compute 2PCFs
    self.__xi_z = []            # 2PCF multipoles to be evaluated
    self.__rmax = 150           # maximum r for the 2PCFs
    self.__dr = 5               # bin size of r for the 2PCFs
    self.__xi_ref = None        # reference 2PCF
    self.__xi_col = []          # columns of the reference 2PCF multipoles

    self.__bk = False           # indicate whether to compute bispectrum
    self.__bkgrid = 256         # grid size for bispectra evaluation
    self.__k1 = [None, None]    # range of the k1 vector for the bispectrum
    self.__k2 = [None, None]    # range of the k2 vector for the bispectrum
    self.__bkbin = 20           # number of output bins for the bispectrum
    self.__bk_ref = None        # reference bispectrum
    self.__bk_col = None        # column of the reference bispectrum

    # Other parameters
    self.__history = []         # history of evaluated parameter sets
    self.__odir = None          # output directory for the current run
    self.__script = None        # script for the current run
    self.__ezfile = None        # EZmock catalogue with RSD for the current run


  def set_param(self, boxsize, num_grid, redshift, num_tracer,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=1, omega_m=0.307115,
      init_pk='/global/u2/z/zhaoc/work/pyEZmock/data/PlanckDM.linear.pk'):
    """
    Set parameters for EZmock evaluation.

    Parameters
    ----------
    boxsize: float
        Side length of the cubic periodic box.
    num_grid: int
        Number of grids per side for the density field.
    redshift: float
        Final redshift of the catalogue.
    num_tracer: int
        Number of tracers to be generated.
    pdf_base: float, optional
        Base number for PDF mapping.
    dens_scat: float, optional
        Density modification parameter.
    rand_motion: float, optional
        Parameter for the random local motion.
    dens_cut: float, optional
        The critical density.
    seed: int, optional
        Random seed.
    omega_m: float, optional
        Density parameter at z = 0.
    init_pk: str, optional
        Initial power spectrum.
    Reference
    ---------
    https://arxiv.org/abs/2007.08997
    """
    self.__param['boxsize'] = float(boxsize)
    self.__param['num_grid'] = int(num_grid)
    self.__param['redshift'] = float(redshift)
    self.__param['num_tracer'] = int(num_tracer)
    if not pdf_base is None: self.__param['pdf_base'] = float(pdf_base)
    if not dens_scat is None: self.__param['dens_scat'] = float(dens_scat)
    if not rand_motion is None: self.__param['rand_motion'] = float(rand_motion)
    if not dens_cut is None: self.__param['dens_cut'] = float(dens_cut)
    self.__param['seed'] = int(seed)

    if not init_pk is None: self.__param['init_pk'] = init_pk
    if not os.path.isfile(self.__param['init_pk']):
      raise IOError(f"init_pk does not exist: {self.__param['init_pk']}")
    self.__param['init_pk'] = os.path.abspath(self.__param['init_pk'])

    self.__param['omega_m'] = float(omega_m)
    if self.__param['omega_m'] <= 0 or self.__param['omega_m'] > 1:
      raise ValueError('omega_m must be between 0 and 1')


  def set_clustering(self, pk=False, pk_ell=[0,2], pk_grid=256, pk_kmax=0.3,
      pk_dk=0.01, pk_ref=None, pk_ref_col=[5,6], xi=False, xi_ell=[0,2],
      xi_rmax=150, xi_dr=5, xi_ref=None, xi_ref_col=[4,5], bk=False,
      bk_grid=256, bk_k1=[0.04,0.06], bk_k2=[0.09,0.11], bk_nbin=20,
      bk_ref=None, bk_ref_col=4):
    """
    Set configurations for clustering measurements.

    Parameters
    ----------
    pk: bool, optional
        Indicate whether to compute power spectrum.
    pk_ell: tuple of ints, optinal
        Legendre multipoles of power spectrum to be evaluated.
    pk_grid: int, optional
        Grid size for power spectra evaluations.
    pk_kmax: float, optional
        Maximum k for the power spectra.
    pk_dk: float, optional
        Bin size of k for the power spectra.
    pk_ref: str, optional
        File for the reference power spectrum. The first column must be k.
    pk_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference power spectrum multipoles.

    xi: bool, optional
        Indicate whether to compute 2-point correlation function (2PCF).
    xi_ell: tuple of ints, optinal
        Legendre multipoles of 2PCF to be evaluated.
    xi_rmax: float, optional
        Maximum separation for the 2PCFs.
    xi_dr: float, optional
        Bin size of separation for the 2PCFs.
    xi_ref: str, optional
        File for the reference 2PCF. The first column must be r.
    xi_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference 2PCF multipoles.

    bk: bool, optional
        Indicate whether to compute bispectrum.
    bk_grid: int, optional
        Grid size for bispectra evaluations.
    bk_k1: (float, float), optional
        Range of the first k vector for the bispectra.
    bk_k2: (float, float), optional
        Range of the second k vector for the bispectra.
    bk_nbin: int, optional
        Number of output bins for the bispectra.
    bk_ref: str, optional
        File for the reference bispectrum. The first column must be theta.
    bk_ref_col: int, optional
        Column of the reference bispectrum.
    """
    self.__pk = pk
    if self.__pk and self.pk_exe is None:
      raise ValueError('pk_exe must be available if computing power spectrum')
    if self.__pk and len(pk_ell) == 0: raise ValueError('pk_ell is empty')
    self.__pkl = pk_ell
    self.__pkgrid = pk_grid
    self.__kmax = pk_kmax
    self.__dk = pk_dk
    if (not pk_ref is None) and (not os.path.isfile(pk_ref)):
      raise IOError(f'pk_ref does not exist: {pk_ref}')
    self.__pk_ref = pk_ref
    if not pk_ref is None and len(pk_ell) != len(pk_ref_col):
      raise ValueError('pk_ell and pk_ref_col should have equal length')
    self.__pk_col = pk_ref_col

    self.__xi = xi
    if self.__xi and self.xi_exe is None:
      raise ValueError('xi_exe must be available if computing 2PCF')
    if self.__xi and len(xi_ell) == 0: raise ValueError('xi_ell is empty')
    self.__xil = xi_ell
    self.__rmax = xi_rmax
    self.__dr = xi_dr
    if (not xi_ref is None) and (not os.path.isfile(xi_ref)):
      raise IOError(f'xi_ref does not exist: {xi_ref}')
    self.__xi_ref = xi_ref
    if not xi_ref is None and len(xi_ell) != len(xi_ref_col):
      raise ValueError('xi_ell and xi_ref_col should have equal length')
    self.__xi_col = xi_ref_col

    self.__bk = bk
    if self.__bk and self.bk_exe is None:
      raise ValueError('bk_exe must be available if computing bispectrum')
    self.__bkgrid = bk_grid
    self.__k1 = bk_k1
    self.__k2 = bk_k2
    self.__bkbin = bk_nbin
    if (not bk_ref is None) and (not os.path.isfile(bk_ref)):
      raise IOError(f'bk_ref does not ebkst: {bk_ref}')
    self.__bk_ref = bk_ref
    self.__bk_col = int(bk_ref_col)


  def run(self, nthreads, queue=None, walltime=30, partition='haswell',
      boxsize=None, num_grid=None, redshift=None, num_tracer=None,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=None, omega_m=None, init_pk=None):
    """
    Run the job for EZmock generation and clustering measurements.

    Parameters
    ----------
    nthreads: int
        Number of OpenMP threads used for the run.
    queue: str, optional
        Queue of the job to be submitted to (e.g. 'debug' and 'regular').
        If not provided, the job script has to be run manually.
    walltime: int, optional
        Limit on the total run time (in minutes) of the job.
    partition: str, optional
        Specify the architecture of the nodes for the job.
        It has to be 'haswell' or 'knl'.
    The rest of the parameters are the same as those for `set_param`.

    Return
    ------
    Filename of the jobs cript.
    """
    import copy
    from subprocess import Popen, PIPE

    nthreads = int(nthreads)
    if nthreads <= 0: raise ValueError(f'Invalid nthreads: {nthreads:d}')
    if not queue is None:
      if walltime is None:
        raise ValueError('walltime is required if queue = None')
      if partition != 'haswell' and partition != 'knl':
        raise ValueError('Invalid job partition: {partition}')
    previous_param = copy.deepcopy(self.__param)

    if not boxsize is None: self.__param['boxsize'] = float(boxsize)
    if not num_grid is None: self.__param['num_grid'] = int(num_grid)
    if not redshift is None: self.__param['redshift'] = float(redshift)
    if not num_tracer is None: self.__param['num_tracer'] = int(num_tracer)
    if not pdf_base is None: self.__param['pdf_base'] = float(pdf_base)
    if not dens_scat is None: self.__param['dens_scat'] = float(dens_scat)
    if not rand_motion is None: self.__param['rand_motion'] = float(rand_motion)
    if not dens_cut is None: self.__param['dens_cut'] = float(dens_cut)
    if not seed is None: self.__param['seed'] = int(seed)
    if not omega_m is None:
      self.__param['omega_m'] = float(omega_m)
      if self.__param['omega_m'] <= 0 or self.__param['omega_m'] > 1:
        raise ValueError('omega_m must be between 0 and 1')
    if not init_pk is None:
      if not os.path.isfile(init_pk):
        raise IOError(f'init_pk does not exist: {init_pk}')
      self.__param['init_pk'] = init_pk

    if None in self.__param.values():
      raise ValueError('Please set EZmock parameters via `set_param` or `run`.')

    # Create the path and name of the job script for the current run
    bname = self.__get_bname(self.__param)
    self.__odir = f'{self.workdir}/{bname}'
    if not os.path.isdir(self.__odir):
      try: os.mkdir(self.__odir)
      except: raise IOError(f'Cannot create directory: {self.__odir}')
    self.__script = f'{self.__odir}/run_{bname}.sh'
    if os.path.isfile(self.__script):
      raise IOError('The job script exists, please wait if the job ' \
          'has already been submitted, or submit/run the script manually: \n'
          f'{self.__script}')
    self.__ezfile = f'EZmock_{bname}_RSD.dat'

    # Store the previous set of parameters as history
    if not None in previous_param.values():
      prev_bname = self.__get_bname(previous_param)
      if os.path.isfile(f'{self.workdir}/{prev_bname}/DONE'):
        self.__history.append(previous_param)

    # Generate contents of the job script file.
    jobstr = ('#!/bin/bash\n#SBATCH -n 1\n#SBATCH -L SCRATCH\n'
        f'#SBATCH -o {self.__odir}/stdout_%j.txt\n'
        f'#SBATCH -e {self.__odir}/stderr_%j.txt\n')
    if not queue is None:
      jobstr += f'#SBATCH -q {queue}\n'
      jobstr += f'#SBATCH -C {partition}\n#SBATCH -c {nthreads:d}\n'
      jobstr += '#SBATCH -t {:d}\n'.format(int(walltime))

    jobstr += f'\nexport OMP_NUM_THREADS={nthreads:d}\n\ncd {self.__odir}\n\n'

    run_mock = True
    if os.path.isfile(f'{self.__odir}/{self.__ezfile}'):
      warn(('EZmock will not be run, as file exists: '
            f'{self.__odir}/{self.__ezfile}'))
      run_mock = False
    if run_mock: jobstr += self.__mock_cmd(bname)

    if self.__pk: jobstr += self.__pk_cmd()
    if self.__xi: jobstr += self.__xi_cmd()
    if self.__bk: jobstr += self.__bk_cmd()

    jobstr += 'echo 1 > DONE\n'

    # Save the job script and submit it if applicable
    with open(self.__script, 'w') as f: f.write(jobstr)

    if queue is None:
      print('Job script generated. Please run the following command manually:')
      print(f' \nbash {self.__script}')
    else:
      process = Popen(['/usr/bin/sbatch',self.__script], shell=False,
          stdout=PIPE, stderr=PIPE, text=True)
      sts = process.wait()
      for line in process.stdout: print(line, end='')
      for line in process.stderr: print(line, end='')
      if sts != 0:
        print('Job submission failed. Please resubmit the script manually:')
        print(f'sbatch {self.__script}')
        print(f'or run the script directly:\nbash {self.__script}')

    return self.__script


  def plot(self, fname=None):
    """
    Plot the clustering measurements of the previous runs, the references,
    and the current run.

    Parameters
    ----------
    fname: str, optional
        If fname is provided, the plot is saved as this file.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Determine the number of subplots
    nplot = 0
    if self.__pk: nplot += len(self.__pkl)
    if self.__xi: nplot += len(self.__xil)
    if self.__bk: nplot += 1
    if nplot == 0:
      print('No clustering measurements specified.')
      print('Please set via `set_clustering`')
      return

    # Create subplots
    ncol = 3
    if nplot <= 4: ncol = nplot
    nrow = int(np.ceil(nplot / ncol))
    figw = min(15, nplot * 5)
    figh = 3 * nrow
    plt.figure(figsize=(figw,figh))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax = [plt.subplot2grid((nrow,ncol), (i//ncol,i%ncol)) for i in range(nplot)]
    for a in ax: a.grid(ls=':', c='dimgray', alpha=0.6)

    iplot = 0
    alpha_hist = 0.8
    ls_ref = 'k:'
    ls_curr = 'k-'

    # Plot power spectra multiples
    if self.__pk:
      # Plot histories first
      for j, hist in enumerate(self.__history):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/PK_EZmock_{bname}_RSD.dat'
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.__pkl):
          ax[iplot+i].plot(d[0], d[5+i]*d[0]**1.5, alpha=alpha_hist,
              label=f'history {j:d}')
      # Plot the reference
      if not self.__pk_ref is None:
        d = np.loadtxt(self.__pk_ref, unpack=True)
        for i, c in enumerate(self.__pk_col):
          ax[iplot+i].plot(d[0], d[c]*d[0]**1.5, ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.__odir}/PK_{self.__ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.__pkl):
          ax[iplot+i].plot(d[0], d[5+i]*d[0]**1.5, ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels and ranges
      for i, ell in enumerate(self.__pkl):
        ax[iplot+i].set_xlabel(r'$k$')
        ax[iplot+i].set_ylabel(r'$k^{{1.5}} P_{:d} (k)$'.format(ell))
        ax[iplot+i].set_xlim(0, self.__kmax)
      iplot += len(self.__pkl)

    # Plot 2PCF
    if self.__xi:
      # Plot histories first
      for j, hist in enumerate(self.__history):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/2PCF_EZmock_{bname}_RSD.dat'
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.__xil):
          ax[iplot+i].plot(d[0], d[3+i]*d[0]**2, alpha=alpha_hist,
              label=f'history {j:d}')
      # Plot the reference
      if not self.__xi_ref is None:
        d = np.loadtxt(self.__xi_ref, unpack=True)
        for i, c in enumerate(self.__xi_col):
          ax[iplot+i].plot(d[0], d[c]*d[0]**2, ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.__odir}/2PCF_{self.__ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.__xil):
          ax[iplot+i].plot(d[0], d[3+i]*d[0]**2, ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels and ranges
      for i, ell in enumerate(self.__xil):
        ax[iplot+i].set_xlabel(r'$s$')
        ax[iplot+i].set_ylabel(r'$s^{{2}} \xi_{:d} (s)$'.format(ell))
        ax[iplot+i].set_xlim(0, self.__rmax)
      iplot += len(self.__xil)

    if self.__bk:
      # Plot histories first
      for j, hist in enumerate(self.__history):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/BK_EZmock_{bname}_RSD.dat'
        d = np.loadtxt(ifile, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], alpha=alpha_hist,
            label=f'history {j:d}')
      # Plot the reference
      if not self.__bk_ref is None:
        d = np.loadtxt(self.__bk_ref, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.__odir}/BK_{self.__ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels
      ax[iplot].set_xlabel(r'$\theta_{12} / \pi$')
      ax[iplot].set_ylabel(r'$B (\theta_{12})$')
      iplot += 1

    ax[-1].legend()
    if not fname is None: plt.savefig(fname)


  def massive_jobs(self, nthreads, seeds, clustering=False, queue=None,
      walltime=30, partition='haswell'):
    """
    Create the job script for massive production of EZmock,
    with the current set of parameters.

    Parameters
    ----------
    nthreads: int
        Number of OpenMP threads used for each realisation.
    seeds: tuple of int
        Random seeds for the massive production.
    clustering: bool, optional
        Indicate whether to perform clustering measurements.
    queue: str, optional
        Queue of the job to be submitted to (e.g. 'regular')
    walltime: int, optional
        Limit on the total run time (in minutes) of the job.
    partition: str, optional
        Specify the architecture of the nodes for the job.
        It has to be 'haswell' or 'knl'.

    Return
    ------
        Filename of the job script.
    """
    import copy

    nthreads = int(nthreads)
    if nthreads <= 0: raise ValueError(f'Invalid nthreads: {nthreads:d}')
    if not queue is None:
      if walltime is None:
        raise ValueError('walltime is required if queue = None')
      if partition != 'haswell' and partition != 'knl':
        raise ValueError('Invalid job partition: {partition}')
    if (len(seeds) <= 1): raise ValueError('seeds must be a list')
    if None in self.__param.values():
      raise ValueError('Please set EZmock parameters via `set_param` or `run`.')

    jobs = []
    # Generate the script for each seed.
    for s in seeds:
      if s <= 0:
        warn(f'Omitting non-positive seed: {s:d}')
        continue
      params = copy.deepcopy(self.__param)
      params['seed'] = s
      bname = self.__get_bname(params)
      odir = f'{self.workdir}/{bname}'
      if not os.path.isdir(odir):
        try: os.mkdir(odir)
        except: raise IOError(f'Cannot create directory: {odir}')
      script = f'{odir}/run_{bname}.sh'
      if os.path.isfile(script):
        warn(f'Omitting existing job for seed {s:d}:\n{script}')
        continue

      jobstr = (f'#!/bin/bash\n\nexport OMP_NUM_THREADS={nthreads:d}\n\n'
          f'cd {odir}\n\n')
      if clustering:
        jobstr += self.__mock_cmd(bname, params=params)
        if self.__pk: jobstr += self.__pk_cmd(bname)
        if self.__xi: jobstr += self.__xi_cmd(bname)
        if self.__bk: jobstr += self.__bk_cmd(bname)
      else:
        jobstr += self.__mock_cmd(bname, rsd=False, params=params)

      # Save the script
      with open(script, 'w') as f: f.write(jobstr)
      jobs.append(script)

    # Generate the job script for all realisations.
    njob = len(jobs)
    if njob < 1: raise ValueError('No valid job is found')
    script = f'{self.workdir}/submit_massive_production.sh'
    jobstr = '#!/bin/bash\n'
    if not queue is None:
      jobstr += (f'#SBATCH -n {njob:d}\n#SBATCH -L SCRATCH\n'
        f'#SBATCH -o {self.workdir}/massive_stdout_%j.txt\n'
        f'#SBATCH -e {self.workdir}/massive_stderr_%j.txt\n')
      jobstr += (f'#SBATCH -q {queue}\n#SBATCH -C {partition}\n'
        f'#SBATCH -c {nthreads:d}\n#SBATCH -t {int(walltime):d}\n')
      for j in jobs:
        jobstr += f'srun -n 1 -c {nthreads:d} --cpu_bind=cores bash {j} &\n'
      jobstr += 'wait\n'

      with open(script, 'w') as f: f.write(jobstr)
      print(f'The job script for {njob:d} realisations have been generated.\n'
          f'Please check it before submission:\n{script}')
    else:
      for j in jobs:
        jobstr += f'bash {j}\n'
      with open(script, 'w') as f: f.write(jobstr)
      print(f'The job script for {njob:d} realisations have been generated.\n'
          f'Please consider running it with `jobfork`:\n{script}')

    return script


  def params(self):
    """
    Print the current set of EZmock parameters.
    """
    for key, value in self.__param.items(): print(key, '=', value)


  def history(self):
    """
    Print the histories of the EZmock parameter sets.
    """
    for i, param in enumerate(self.__history): print(f'{i:d}:', param)


  def clear(self, slicer):
    """
    Clear history entries defined by `slicer`

    Parameters
    ----------
    slicer:
        The slice of the histories to be cleared.
        It must be generated using the `slice` function.
    """
    if type(slicer).__name__ != 'slice':
      raise TypeError('slicer must be generated using the `slice` function')
    del(self.__history[slicer])


  def __get_bname(self, param):
    """
    Generate the basename of files for a given parameter set.

    Parameters
    ----------
    param: dict
        Dictionary storing a set of EZmock parameters.

    Return
    ------
    The basename as a string.
    """
    if None in param.values(): return None
    bname = (f"B{param['boxsize']:g}"
        f"G{param['num_grid']:d}"
        f"Z{param['redshift']:g}"
        f"N{param['num_tracer']:d}_"
        f"b{param['pdf_base']:g}"
        f"d{param['dens_scat']:g}"
        f"r{param['rand_motion']:g}"
        f"c{param['dens_cut']:g}_"
        f"seed{param['seed']:d}")
    return bname


  def __mock_cmd(self, bname, rsd=True, params=None):
    """
    Generate the command for generating the redshift space EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.
    rsd: bool, optional
        Indicate whether to apply redshift space distortions.
    params: dict, optional
        A dictionary of the EZmock parameters.

    Returns
    -------
    The command as a string
    """
    from cosmoEZ import flatLCDM

    if params is None: params = self.__param

    # Compute structure growth parameters
    cosmo = flatLCDM(omega_m = params['omega_m'])
    z1 = 1 + params['redshift']
    a = 1. / z1
    grow2z0 = cosmo.growth2z0(a)
    hubble = cosmo.hubble(a)
    zdist = cosmo.growthf(a) * hubble * a

    # Generate the command for running EZmock
    jobstr = ("echo '&EZmock_v0_input\n"
        f"datafile_prefix = \"EZmock_{bname}\"\n"
        f"datafile_path = \"./\"\n"
        f"iseed = {params['seed']:d}\n"
        f"boxsize = {params['boxsize']:g}\n"
        f"grid_num = {params['num_grid']:d}\n"
        f"redshift = {params['redshift']:g}\n"
        f"grow2z0 = {grow2z0:g}\n"
        f"expect_sum_pdf = {params['num_tracer']:d}\n"
        f"expect_A_pdf = {params['pdf_base']:g}\n"
        f"density_cut = {params['dens_cut']:g}\n"
        f"scatter2 = {params['dens_scat']:g}\n"
        f"zdist_rate = {zdist:g}\n"
        f"zdist_fog = {params['rand_motion']:g}\n"
        "density_sat = 100\nscatter = 10\nmodify_pk = 0.0\n"
        "modify_pdf = 0\nantidamping = 2\n"
        "use_whitenoise_file = .false.\nwhitenoise_file = \"\"\n"
        f"pkfile = \"{params['init_pk']}\"\n"
        f"pknwfile = \"{params['init_pk']}\"\n"
        "compute_CF = .false.\ncompute_CF_zdist = .false.\n"
        "dilute_factor = 0.3\nskiplines = 0\ntwod_corr_suffix = \"\"\n"
        f"max_r = 50\nbin_size = 5\nom = {params['omega_m']:g}\n/'"
        f" | {self.ez_exe} || exit\n\n")

    if rsd:
      # Generate the command for applying redshift space distortions
      bsize = params['boxsize']
      # Check boundaries and remove non-numerical entries (such as nan)
      jobstr += ("awk '{CONVFMT=\"%.8g\";OFMT=\"%.8g\"; "
          "if ($3+0==$3 && $6+0==$6) { "
          f"z=($3+$6*{z1:g}/{hubble:.8g}+{bsize:g})%{bsize:g}; "
          f"if ($1>=0 && $1<{bsize:g} && $2>=0 && $2<{bsize:g}"
          f" && z>=0 && z<{bsize:g}) print $1,$2,z; }} }} ' "
          f"EZmock_{bname}.dat > EZmock_{bname}_RSD.dat || exit\n\n")

    return jobstr


  def __pk_cmd(self, bname=None):
    """
    Generate the command for computing power spectrum of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.

    Returns
    -------
    The command as a string
    """
    bsize = self.__param['boxsize']
    poles = '[' + ','.join(f'{i:d}' for i in self.__pkl) + ']'
    ifile = self.__ezfile
    if not bname is None: ifile = f'EZmock_{bname}_RSD.dat'
    ofile = f'PK_{ifile}'

    jobstr = (f'{self.pk_exe} -d {ifile} --data-formatter "%lf %lf %lf" '
        f"-p '[($1+{bsize:g})%{bsize:g},($2+{bsize:g})%{bsize:g},"
        f"($3+{bsize:g})%{bsize:g}]' -s T -B {bsize:g} -G {self.__pkgrid:d} "
        f"-n 1 -i F -l '{poles}' -k 0 -K {self.__kmax:g} -b {self.__dk:g} "
        f"-a {ofile} || exit\n\n")

    return jobstr


  def __xi_cmd(self, bname=None):
    """
    Generate the command for computing 2-point correlation function of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.

    Returns
    -------
    The command as a string
    """
    bsize = self.__param['boxsize']
    poles = '[' + ','.join(f'{i:d}' for i in self.__xil) + ']'
    ifile = self.__ezfile
    if not bname is None: ifile = f'EZmock_{bname}_RSD.dat'
    ofile = f'2PCF_{ifile}'

    jobstr = (f"{self.xi_exe} -i {ifile} -l D -f '%lf %lf %lf' -x '[$1,$2,$3]' "
        f'-b {bsize:g} -B 1 -p DD -P {ofile}.dd -e "DD/@@-1" -E {ofile}.xi2d '
        f"-m '{poles}' -M {ofile} --s-min 0 --s-max {self.__rmax:g} "
        f'--s-step {self.__dr:g} --mu-num 60 --dist-prec 0 -S 0 || exit\n\n')

    return jobstr


  def __bk_cmd(self, bname=None):
    """
    Generate the command for computing power spectrum of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.

    Returns
    -------
    The command as a string
    """
    bsize = self.__param['boxsize']
    ifile = self.__ezfile
    if not bname is None: ifile = f'EZmock_{bname}_RSD.dat'
    ofile = f'BK_{ifile}'

    jobstr = (f'{self.bk_exe} -i {ifile} -s 7 --x-min 0 --y-min 0 --z-min 0 '
        f'--x-max {bsize:g} --y-max {bsize:g} --z-max {bsize:g} -b 0'
        f' -B {bsize:g} -g {self.__bkgrid:d} -w 1 -x 0 -p {self.__k1[0]:g} '
        f'-P {self.__k1[1]:g} -q {self.__k2[0]:g} -Q {self.__k2[1]:g} '
        f'-n {self.__bkbin:d} -o {ofile} -y || exit\n\n')

    return jobstr

