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
import json

class pyEZmock:
  """
  Python wrapper of Effective Zel'dovich approximation mock (EZmock).
  """

  def __init__(self, workdir, restore=False,
      exe='/global/u2/z/zhaoc/work/pyEZmock/bin/EZmock.sh',
      pk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/POWSPEC.sh',
      xi_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/FCFC_2PT_BOX.sh',
      bk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/BISPEC_BOX.sh'):
    """
    Initialise the `pyEZmock` class.

    Parameters
    ----------
    workdir: str
        Working directory for generating EZmocks.
    restore: bool, optional
        Indicate whether to scan `workdir` and restore previous runs.
    exe: str, optional
        Location of the EZmock executable.
    pk_exe: str, optional
        Location of the powspec executable.
    xi_exe: str, optional
        Location of the FCFC executable for periodic boxes.
    bk_exe: str, optional
        Location of the bispec executable.
    """
    # Set paths and filenames
    self.workdir = workdir
    self.exe = exe
    self.pk_exe = pk_exe
    self.xi_exe = xi_exe
    self.bk_exe = bk_exe

    self._conf = 'config.json'  # file for storing configurations
    self._script = 'run_job.sh' # job script file
    self._done = 'DONE'         # indicate whether the job is finished

    # Initialise EZmock parameters
    self._param = dict(
      boxsize = None,
      num_grid = None,
      redshift = None,
      num_tracer = None,
      pdf_base = None,
      dens_scat = None,
      rand_motion = None,
      dens_cut = None,
      seed = 1,
      omega_m = None,
      z_init = 0,
      init_pk = None
    )

    # Initialise clustering settings
    self._pkconf = dict(        # Power spectrum (pk) settings
      rspace = False,           # indicate whether to compute real-space pk
      zspace = False,           # indicate whether to compute redshift-space pk
      rell = [],                # real-space pk multipoles to be evaluated
      zell = [],                # redshift-space pk multipoles to be evaluated
      ngrid = 256,              # grid size for power spectra evaluation
      kmax = 0.3,               # maximum k for the power spectra
      dk = 0.01,                # bin size of k for the power spectra
      rref = None,              # reference real-space power spectrum
      rcol = [],                # reference real-space pk multipole columns
      zref = None,              # reference redshift-space power spectrum
      zcol = []                 # reference redshift-space pk multipole columns
    )

    self._xiconf = dict(        # 2-point correlation function (xi) settings
      rspace = False,           # indicate whether to compute real-space xi
      zspace = False,           # indicate whether to compute redshift-space xi
      rell = [],                # real-space xi multipoles to be evaluated
      zell = [],                # redshift-space xi multipoles to be evaluated
      rmax = 150,               # maximum separation for the 2PCFs
      dr = 5,                   # separation bin size for the 2PCFs
      nmu = 60,                 # number of mu bins
      rref = None,              # reference real-space xi
      rcol = [],                # reference real-space xi multipole columns
      zref = None,              # reference redshift-space xi
      zcol = []                 # reference redshift-space xi multipole columns
    )

    self._bkconf = dict(        # Bispectrum (bk) settings
      rspace = False,           # indicate whether to compute real-space bk
      zspace = False,           # indicate whether to compute redshift-space bk
      ngrid = 256,              # grid size for bispectra evaluation
      k1 = [None, None],        # range of the k1 vector for the bispectrum
      k2 = [None, None],        # range of the k2 vector for the bispectrum
      nbin = 20,                # number of output bins for the bispectrum
      rref = None,              # reference real-space bispectrum
      rcol = None,              # column of the reference real-space bk
      zref = None,              # reference redshift-space bispectrum
      zcol = None               # column of the reference redshift-space bk
    )

    # Default columns for clustering measurements
    self._pkcol = 5
    self._xicol = 3
    self._bkcol = 4

    # Default plotting styles
    self._alpha = 0.8           # transparency for historical curves
    self._ls_ref = 'k:'         # line style for the reference
    self._ls_curr = 'k-'        # line style for the current run

    # Other parameters
    self._history = []         # history of evaluated parameter sets
    self._odir = None          # output directory for the current run
    self._bname = None         # catalogue basename for the current run

    # Load previous runs as histories
    if restore: self.restore()


  def set_param(self, boxsize, num_grid, redshift, num_tracer,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=1, omega_m=0.307115, z_init=0,
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
    z_init: float, optional
        Redshift of the initial power spectrum.
    init_pk: str, optional
        Initial power spectrum.
    Reference
    ---------
    https://arxiv.org/abs/2007.08997
    """
    self._param['boxsize'] = float(boxsize)
    self._param['num_grid'] = int(num_grid)
    self._param['redshift'] = float(redshift)
    self._param['num_tracer'] = int(num_tracer)
    if pdf_base is not None: self._param['pdf_base'] = float(pdf_base)
    if dens_scat is not None: self._param['dens_scat'] = float(dens_scat)
    if rand_motion is not None: self._param['rand_motion'] = float(rand_motion)
    if dens_cut is not None: self._param['dens_cut'] = float(dens_cut)
    self._param['seed'] = int(seed)

    if init_pk is not None:
      if not os.path.isfile(init_pk):
        raise FileNotFoundError(f"`init_pk` does not exist: {init_pk}")
      self._param['init_pk'] = os.path.abspath(init_pk)

    self._param['omega_m'] = float(omega_m)
    if self._param['omega_m'] <= 0 or self._param['omega_m'] > 1:
      raise ValueError('`omega_m` must be between 0 and 1')
    self._param['z_init'] = float(z_init)
    if self._param['z_init'] < 0:
      raise ValueError('`z_init` must be non-negative')


  def set_clustering(self,
      pk='none', pk_r_ell=[0], pk_z_ell=[0,2], pk_grid=256, pk_kmax=0.3,
      pk_dk=0.01, pk_r_ref=None, pk_r_ref_col=[5], pk_z_ref=None,
      pk_z_ref_col=[5,6],
      xi='none', xi_r_ell=[0], xi_z_ell=[0,2], xi_rmax=150, xi_dr=5, xi_nmu=60,
      xi_r_ref=None, xi_r_ref_col=[4], xi_z_ref=None, xi_z_ref_col=[4,5],
      bk='none', bk_grid=256, bk_k1=[0.04,0.06], bk_k2=[0.09,0.11], bk_nbin=20,
      bk_r_ref=None, bk_r_ref_col=4, bk_z_ref=None, bk_z_ref_col=4):
    """
    Set configurations for clustering measurements.

    Parameters
    ----------
    pk: str, optional
        Specify the power spectra to be computed. The string has to be one of
        'real', 'redshift', 'both', or 'none'. 'real' and 'redshift' indicate
        computing real- and redshift-space power spectrum multipoles
        respectively; 'both' and 'none' indicate computing both or none of them.
    pk_r_ell: tuple of ints, optinal
        Legendre multipoles of real-space power spectrum to be evaluated.
    pk_z_ell: tuple of ints, optinal
        Legendre multipoles of redshift-space power spectrum to be evaluated.
    pk_grid: int, optional
        Grid size for power spectra evaluations.
    pk_kmax: float, optional
        Maximum k for the power spectra.
    pk_dk: float, optional
        Bin size of k for the power spectra.
    pk_r_ref: str, optional
        File for the reference real-space power spectrum.
        The first column must be k.
    pk_r_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference real-space power spectrum
        multipoles.
    pk_z_ref: str, optional
        File for the reference redshift-space power spectrum.
        The first column must be k.
    pk_z_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference redshift-space power spectrum
        multipoles.

    xi: str, optional
        Specify the 2-point correlation functions (2PCFs) to be computed.
        The string has to be one of 'real', 'redshift', 'both', or 'none'.
    xi_r_ell: tuple of ints, optinal
        Legendre multipoles of real-space 2PCF to be evaluated.
    xi_z_ell: tuple of ints, optinal
        Legendre multipoles of redshift-space 2PCF to be evaluated.
    xi_rmax: float, optional
        Maximum separation for the 2PCFs.
    xi_dr: float, optional
        Bin size of separation for the 2PCFs.
    xi_nmu: int, optional
        Number of mu bins for the 2PCFs.
    xi_r_ref: str, optional
        File for the reference real-space 2PCF. The first column must be r.
    xi_r_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference real-space 2PCF multipoles.
    xi_z_ref: str, optional
        File for the reference redshift-space 2PCF. The first column must be r.
    xi_z_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference redshift-space 2PCF
        multipoles.

    bk: str, optional
        Specify the bispectra to be computed. The string has to be one of
        'real', 'redshift', 'both', or 'none'.
    bk_grid: int, optional
        Grid size for bispectra evaluations.
    bk_k1: (float, float), optional
        Range of the first k vector for the bispectra.
    bk_k2: (float, float), optional
        Range of the second k vector for the bispectra.
    bk_nbin: int, optional
        Number of output bins for the bispectra.
    bk_r_ref: str, optional
        File for the reference real-space bispectrum.
        The first column must be theta.
    bk_r_ref_col: int, optional
        Column (counting from 0) of the reference real-space bispectrum.
    bk_z_ref: str, optional
        File for the reference redshift-space bispectrum.
        The first column must be theta.
    bk_z_ref_col: int, optional
        Column (counting from 0) of the reference redshift-space bispectrum.
    """
    if self._pk_exe is not None:
      if pk not in ('real', 'redshift', 'both', 'none'):
        raise ValueError("`pk` must be 'real', 'redshift', 'both' or 'none'")
      self._pkconf['rspace'] = (pk == 'real') or (pk == 'both')
      self._pkconf['zspace'] = (pk == 'redshift') or (pk == 'both')
      if pk != 'none':
        self._pkconf['ngrid'] = int(pk_grid)
        self._pkconf['kmax'] = float(pk_kmax)
        self._pkconf['dk'] = float(pk_dk)
      if self._pkconf['rspace']:
        try:
          if len(pk_r_ell) == 0: raise ValueError('`pk_r_ell` is empty')
        except TypeError: pk_r_ell = [pk_r_ell]
        self._pkconf['rell'] = pk_r_ell
        if pk_r_ref is not None:
          if not os.path.isfile(pk_r_ref):
            raise FileNotFoundError(f'`pk_r_ref` does not exist: {pk_r_ref}')
          if len(pk_r_ell) != len(pk_r_ref_col):
            raise ValueError('unequal length of `pk_r_ell` and `pk_r_ref_col`')
          self._pkconf['rref'] = os.path.abspath(pk_r_ref)
          self._pkconf['rcol'] = pk_r_ref_col
      if self._pkconf['zspace']:
        try:
          if len(pk_z_ell) == 0: raise ValueError('`pk_z_ell` is empty')
        except TypeError: pk_z_ell = [pk_z_ell]
        self._pkconf['zell'] = pk_z_ell
        if pk_z_ref is not None:
          if not os.path.isfile(pk_z_ref):
            raise FileNotFoundError(f'`pk_z_ref` does not exist: {pk_z_ref}')
          if len(pk_z_ell) != len(pk_z_ref_col):
            raise ValueError('unequal length of `pk_z_ell` and `pk_z_ref_col`')
          self._pkconf['zref'] = os.path.abspath(pk_z_ref)
          self._pkconf['zcol'] = pk_z_ref_col
    else:
      self._pkconf['rspace'] = self._pkconf['zspace'] = False

    if self._xi_exe is not None:
      if xi not in ('real', 'redshift', 'both', 'none'):
        raise ValueError("`xi` must be 'real', 'redshift', 'both' or 'none'")
      self._xiconf['rspace'] = (xi == 'real') or (xi == 'both')
      self._xiconf['zspace'] = (xi == 'redshift') or (xi == 'both')
      if xi != 'none':
        self._xiconf['rmax'] = float(xi_rmax)
        self._xiconf['dr'] = float(xi_dr)
        self._xiconf['nmu'] = int(xi_nmu)
      if self._xiconf['rspace']:
        try:
          if len(xi_r_ell) == 0: raise ValueError('`xi_r_ell` is empty')
        except TypeError: xi_r_ell = [xi_r_ell]
        self._xiconf['rell'] = xi_r_ell
        if xi_r_ref is not None:
          if not os.path.isfile(xi_r_ref):
            raise FileNotFoundError(f'`xi_r_ref` does not exist: {xi_r_ref}')
          if len(xi_r_ell) != len(xi_r_ref_col):
            raise ValueError('unequal length of `xi_r_ell` and `xi_r_ref_col`')
          self._xiconf['rref'] = os.path.abspath(xi_r_ref)
          self._xiconf['rcol'] = xi_r_ref_col
      if self._xiconf['zspace']:
        try:
          if len(xi_z_ell) == 0: raise ValueError('`xi_z_ell` is empty')
        except TypeError: xi_z_ell = [xi_z_ell]
        self._xiconf['zell'] = xi_z_ell
        if xi_z_ref is not None:
          if not os.path.isfile(xi_z_ref):
            raise FileNotFoundError(f'`xi_z_ref` does not exist: {xi_z_ref}')
          if len(xi_z_ell) != len(xi_z_ref_col):
            raise ValueError('unequal length of `xi_z_ell` and `xi_z_ref_col`')
          self._xiconf['zref'] = os.path.abspath(xi_z_ref)
          self._xiconf['zcol'] = xi_z_ref_col
    else:
      self._xiconf['rspace'] = self._xiconf['zspace'] = False

    if self._bk_exe is not None:
      if bk not in ('real', 'redshift', 'both', 'none'):
        raise ValueError("`bk` must be 'real', 'redshift', 'both' or 'none'")
      self._bkconf['rspace'] = (bk == 'real') or (bk == 'both')
      self._bkconf['zspace'] = (bk == 'redshift') or (bk == 'both')
      if bk != 'none':
        self._bkconf['ngrid'] = int(bk_grid)
        if len(bk_k1) != 2: raise ValueError('invalid length of `bk_k1`')
        if len(bk_k2) != 2: raise ValueError('invalid length of `bk_k2`')
        self._bkconf['k1'] = bk_k1
        self._bkconf['k2'] = bk_k2
        self._bkconf['nbin'] = int(bk_nbin)
      if self._bkconf['rspace']:
        if bk_r_ref is not None:
          if not os.path.isfile(bk_r_ref):
            raise FileNotFoundError(f'`bk_r_ref` does not exist: {bk_r_ref}')
          self._bkconf['rref'] = os.path.abspath(bk_r_ref)
          self._bkconf['rcol'] = int(bk_r_ref_col)
      if self._bkconf['zspace']:
        if bk_z_ref is not None:
          if not os.path.isfile(bk_z_ref):
            raise FileNotFoundError(f'`bk_z_ref` does not exist: {bk_z_ref}')
          self._bkconf['zref'] = os.path.abspath(bk_z_ref)
          self._bkconf['zcol'] = int(bk_z_ref_col)
    else:
      self._bkconf['rspace'] = self._bkconf['zspace'] = False


  def run(self, nthreads, queue=None, walltime=30, partition='haswell',
      boxsize=None, num_grid=None, redshift=None, num_tracer=None,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=None, omega_m=None, z_init=None, init_pk=None):
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
    if nthreads <= 0: raise ValueError(f'invalid `nthreads`: {nthreads:d}')
    if queue is not None:
      if walltime is None:
        raise ValueError('`walltime` is required when `queue` is set')
      if partition != 'haswell' and partition != 'knl':
        raise ValueError("`partition` must be 'haswell' or 'knl'")
    previous_param = copy.deepcopy(self._param)

    if boxsize is not None: self._param['boxsize'] = float(boxsize)
    if num_grid is not None: self._param['num_grid'] = int(num_grid)
    if redshift is not None: self._param['redshift'] = float(redshift)
    if num_tracer is not None: self._param['num_tracer'] = int(num_tracer)
    if pdf_base is not None: self._param['pdf_base'] = float(pdf_base)
    if dens_scat is not None: self._param['dens_scat'] = float(dens_scat)
    if rand_motion is not None: self._param['rand_motion'] = float(rand_motion)
    if dens_cut is not None: self._param['dens_cut'] = float(dens_cut)
    if seed is not None: self._param['seed'] = int(seed)
    if omega_m is not None:
      self._param['omega_m'] = float(omega_m)
      if self._param['omega_m'] <= 0 or self._param['omega_m'] > 1:
        raise ValueError('`omega_m` must be between 0 and 1')
    if z_init is not None:
      self._param['z_init'] = float(z_init)
      if self._param['z_init'] < 0:
        raise ValueError('`z_init` must be non-negative')
    if init_pk is not None:
      if not os.path.isfile(init_pk):
        raise FileNotFoundError(f'`init_pk` does not exist: {init_pk}')
      self._param['init_pk'] = os.path.abspath(init_pk)

    if None in self._param.values():
      raise ValueError('please set EZmock parameters via `set_param` or `run`')

    # Create the path for the current run
    self._bname = self._get_bname(self._param)
    self._odir = f'{self._workdir}/{self._bname}'
    if not os.path.isdir(self._odir):
      try: os.mkdir(self._odir)
      except: raise IOError(f'cannot create directory: {self._odir}')

    # Check if the job exists already
    quit = False
    par = pkconf = xiconf = bkconf = None
    script = f'{self._odir}/{self._script}'
    conf = f'{self._odir}/{self._conf}'
    done = f'{self._odir}/{self._done}'
    if os.path.isfile(script):
      try:
        with open(conf, 'r') as f:
          par, pkconf, xiconf, bkconf = json.load(f)
        if par == self._param and pkconf == self._pkconf and \
            xiconf == self._xiconf and bkconf == self._bkconf:
          if os.path.isfile(done): print('The job has been finished already')
          else: self._warn(('the job exists but have not been finished, please '
                  'wait if the job has already been submitted, or submit/run '
                  f'the script manually: \n{script}'))
          quit = True
        elif par is not None and par != self._param:
          os.rename(conf, f'{conf.old}')
          self._warn(('existing EZmock will be overwritten, due to the '
                'change of parameters. The previous settings are moved to \n'
                f'{conf}.old'))
      except FileNotFoundError:
        self._warn(('existing EZmock run detected but the settings are '
              'not found, rerun anyway'))
      if quit == False and os.path.isfile(done): os.remove(done)

    # Check the previous set of parameters, and record as history if applicable
    if not None in previous_param.values():
      prev_bname = self._get_bname(previous_param)
      if not os.path.isfile(f'{self._workdir}/{prev_bname}/{self._done}'):
        self._warn('the previous run may have not been finished')
      elif previous_param not in self._history:
        self._history.append(previous_param)
    if quit: return script

    # Generate contents of the job script file.
    jobstr = ('#!/bin/bash\n#SBATCH -n 1\n#SBATCH -L SCRATCH\n'
        f'#SBATCH -o {self._odir}/stdout_%j.txt\n'
        f'#SBATCH -e {self._odir}/stderr_%j.txt\n')
    if queue is not None:
      jobstr += (f'#SBATCH -q {queue}\n'
          f'#SBATCH -C {partition}\n#SBATCH -c {nthreads:d}\n'
          f'#SBATCH -t {int(walltime):d}\n')

    jobstr += f'\nexport OMP_NUM_THREADS={nthreads:d}\n\ncd {self._odir}\n\n'

    run_mock = True
    ofile = f'{self._odir}/EZmock_{self._bname}.dat'
    if par == self._param and os.path.isfile(ofile):
      self._warn('EZmock will not be run, as file exists: {ofile}')
      run_mock = False
    else:       # remeasure clustering measurements if rerunning EZmock
      pkconf = xiconf = bkconf = None

    run_rsd = self._pkconf['zspace'] or self._xiconf['zspace'] or \
              self._bkconf['zspace']
    ofile = f'{self._odir}/EZmock_{self._bname}_RSD.dat'
    if run_rsd and par == self._param and os.path.isfile(ofile): run_rsd = False
    if run_mock or run_rsd:
      jobstr += self._mock_cmd(self._bname, mock=run_mock, rsd=run_rsd)

    if self._pkconf['rspace']:
      ofile = f'{self._odir}/PK_EZmock_{self._bname}.dat'
      if pkconf != self._pkconf or not os.path.isfile(ofile):
        jobstr += self._pk_cmd(rsd=False)
    if self._pkconf['zspace']:
      ofile = f'{self._odir}/PK_EZmock_{self._bname}_RSD.dat'
      if pkconf != self._pkconf or not os.path.isfile(ofile):
        jobstr += self._pk_cmd(rsd=True)
    if self._xiconf['rspace']:
      ofile = f'{self._odir}/2PCF_EZmock_{self._bname}.dat'
      if xiconf != self._xiconf or not os.path.isfile(ofile):
        jobstr += self._xi_cmd(rsd=False)
    if self._xiconf['zspace']:
      ofile = f'{self._odir}/2PCF_EZmock_{self._bname}_RSD.dat'
      if xiconf != self._xiconf or not os.path.isfile(ofile):
        jobstr += self._xi_cmd(rsd=True)
    if self._bkconf['rspace']:
      ofile = f'{self._odir}/BK_EZmock_{self._bname}.dat'
      if bkconf != self._bkconf or not os.path.isfile(ofile):
        jobstr += self._bk_cmd(rsd=False)
    if self._bkconf['zspace']:
      ofile = f'{self._odir}/BK_EZmock_{self._bname}_RSD.dat'
      if bkconf != self._bkconf or not os.path.isfile(ofile):
        jobstr += self._bk_cmd(rsd=True)

    jobstr += f'echo 1 > {self._done}\n'

    # Save the job script and configurations
    with open(script, 'w') as f: f.write(jobstr)
    with open(conf, 'w') as f:
      json.dump([self._param, self._pkconf, self._xiconf, self._bkconf], f,
          indent=2)

    # Submit the job if applicable
    if queue is None:
      print(('Job script generated. Please run the following command manually:'
            f'\nbash {script}'))
    else:
      process = Popen(['/usr/bin/sbatch',script], shell=False,
          stdout=PIPE, stderr=PIPE, text=True)
      sts = process.wait()
      for line in process.stdout: print(line, end='')
      for line in process.stderr: print(line, end='')
      if sts != 0:
        self._warn(('job submission failed. Please resubmit the script '
              f'manually: \nsbatch {script} \nor run the script directly: \n'
              f'bash {script}'))

    return script


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
    if self._pkconf['rspace']: nplot += len(self._pkconf['rell'])
    if self._pkconf['zspace']: nplot += len(self._pkconf['zell'])
    if self._xiconf['rspace']: nplot += len(self._xiconf['rell'])
    if self._xiconf['zspace']: nplot += len(self._xiconf['zell'])
    if self._bkconf['rspace']: nplot += 1
    if self._bkconf['zspace']: nplot += 1
    if nplot == 0:
      raise ValueError(('no clustering measurements specified, '
            'please set via `set_clustering`'))

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
    if self._pkconf['rspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/PK_EZmock_{bname}.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._pkconf['rell']):
            ax[iplot+i].plot(d[0], d[self._pkcol+i]*d[0]**1.5,
                alpha=self._alpha, label=f'history {j:d}')
      # Plot the reference
      if self._pkconf['rref'] is not None:
        d = np.loadtxt(self._pkconf['rref'], unpack=True)
        for i, c in enumerate(self._pkconf['rcol']):
          ax[iplot+i].plot(d[0], d[c]*d[0]**1.5, self._ls_ref, label='ref')
      # Plot the current results
      if not None in self._param.values():
        ifile = f'{self._workdir}/{self._bname}/PK_EZmock_{self._bname}.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._pkconf['rell']):
            ax[iplot+i].plot(d[0], d[self._pkcol+i]*d[0]**1.5,
                self._ls_curr, label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels and ranges
      for i, ell in enumerate(self._pkconf['rell']):
        ax[iplot+i].set_xlabel(r'$k$ (real space)')
        ax[iplot+i].set_ylabel(r'$k^{{1.5}} P_{:d} (k)$'.format(ell))
        ax[iplot+i].set_xlim(0, self._pkconf['kmax'])
      iplot += len(self._pkconf['rell'])

    if self._pkconf['zspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/PK_EZmock_{bname}_RSD.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._pkconf['zell']):
            ax[iplot+i].plot(d[0], d[self._pkcol+i]*d[0]**1.5,
                alpha=self._alpha, label=f'history {j:d}')
      # Plot the reference
      if self._pkconf['zref'] is not None:
        d = np.loadtxt(self._pkconf['zref'], unpack=True)
        for i, c in enumerate(self._pkconf['zcol']):
          ax[iplot+i].plot(d[0], d[c]*d[0]**1.5, self._ls_ref, label='ref')
      # Plot the current results
      if not None in self._param.values():
        ifile = f'{self._workdir}/{self._bname}/PK_EZmock_{self._bname}_RSD.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._pkconf['zell']):
            ax[iplot+i].plot(d[0], d[self._pkcol+i]*d[0]**1.5,
                self._ls_curr, label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels and ranges
      for i, ell in enumerate(self._pkconf['zell']):
        ax[iplot+i].set_xlabel(r'$k$ (redshift space)')
        ax[iplot+i].set_ylabel(r'$k^{{1.5}} P_{:d} (k)$'.format(ell))
        ax[iplot+i].set_xlim(0, self._pkconf['kmax'])
      iplot += len(self._pkconf['zell'])

    # Plot 2PCF multipoles
    if self._xiconf['rspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/2PCF_EZmock_{bname}.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._xiconf['rell']):
            ax[iplot+i].plot(d[0], d[self._xicol+i]*d[0]**2,
                alpha=self._alpha, label=f'history {j:d}')
      # Plot the reference
      if self._xiconf['rref'] is not None:
        d = np.loadtxt(self._xiconf['rref'], unpack=True)
        for i, c in enumerate(self._xiconf['rcol']):
          ax[iplot+i].plot(d[0], d[c]*d[0]**2, self._ls_ref, label='ref')
      # Plot the current results
      if not None in self._param.values():
        ifile = f'{self._workdir}/{self._bname}/2PCF_EZmock_{self._bname}.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._xiconf['rell']):
            ax[iplot+i].plot(d[0], d[self._xicol+i]*d[0]**2,
                self._ls_curr, label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels and ranges
      for i, ell in enumerate(self._xiconf['rell']):
        ax[iplot+i].set_xlabel(r'$r$ (real space)')
        ax[iplot+i].set_ylabel(r'$r^{{2}} \xi_{:d} (r)$'.format(ell))
        ax[iplot+i].set_xlim(0, self._xiconf['rmax'])
      iplot += len(self._xiconf['rell'])

    if self._xiconf['zspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/2PCF_EZmock_{bname}_RSD.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._xiconf['zell']):
            ax[iplot+i].plot(d[0], d[self._xicol+i]*d[0]**2,
                alpha=self._alpha, label=f'history {j:d}')
      # Plot the reference
      if self._xiconf['zref'] is not None:
        d = np.loadtxt(self._xiconf['zref'], unpack=True)
        for i, c in enumerate(self._xiconf['zcol']):
          ax[iplot+i].plot(d[0], d[c]*d[0]**2, self._ls_ref, label='ref')
      # Plot the current results
      if not None in self._param.values():
        bname = self._bname
        ifile = f'{self._workdir}/{bname}/2PCF_EZmock_{bname}_RSD.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          for i, ell in enumerate(self._xiconf['zell']):
            ax[iplot+i].plot(d[0], d[self._xicol+i]*d[0]**2,
                self._ls_curr, label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels and ranges
      for i, ell in enumerate(self._xiconf['zell']):
        ax[iplot+i].set_xlabel(r'$s$ (redshift space)')
        ax[iplot+i].set_ylabel(r'$s^{{2}} \xi_{:d} (s)$'.format(ell))
        ax[iplot+i].set_xlim(0, self._xiconf['rmax'])
      iplot += len(self._xiconf['zell'])

    # Plot bispectra
    if self._bkconf['rspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/BK_EZmock_{bname}.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          ax[iplot].plot(d[0]/np.pi, d[self._bkcol], alpha=self._alpha,
              label=f'history {j:d}')
      # Plot the reference
      if self._bkconf['rref'] is not None:
        d = np.loadtxt(self._bkconf['rref'], unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[self._bkconf['rcol']], self._ls_ref,
            label='ref')
      # Plot the current results
      if not None in self._param.values():
        ifile = f'{self._workdir}/{self._bname}/BK_EZmock_{self._bname}.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          ax[iplot].plot(d[0]/np.pi, d[self._bkcol], self._ls_curr,
              label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels
      ax[iplot].set_xlabel(r'$\theta_{12} / \pi$ (real space)')
      ax[iplot].set_ylabel(r'$B (\theta_{12})$')
      iplot += 1

    if self._bkconf['zspace']:
      # Plot histories first
      for j, hist in enumerate(self._history):
        bname = self._get_bname(hist)
        ifile = f'{self._workdir}/{bname}/BK_EZmock_{bname}_RSD.dat'
        if os.path.isfile(ifile):
          d = np.loadtxt(ifile, unpack=True)
          ax[iplot].plot(d[0]/np.pi, d[self._bkcol], alpha=self._alpha,
              label=f'history {j:d}')
      # Plot the reference
      if self._bkconf['zref'] is not None:
        d = np.loadtxt(self._bkconf['zref'], unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[self._bkconf['zcol']], self._ls_ref,
            label='ref')
      # Plot the current results
      if not None in self._param.values():
        ifile = f'{self._workdir}/{self._bname}/BK_EZmock_{self._bname}_RSD.dat'
        try:
          d = np.loadtxt(ifile, unpack=True)
          ax[iplot].plot(d[0]/np.pi, d[self._bkcol], self._ls_curr,
              label='current')
        except FileNotFoundError:
          self._warn('the current job may have not been finished')
      else: self._warn('the current job may have not been initialised')
      # Set axis labels
      ax[iplot].set_xlabel(r'$\theta_{12} / \pi$ (redshift space)')
      ax[iplot].set_ylabel(r'$B (\theta_{12})$')
      iplot += 1

    ax[-1].legend()
    if fname is not None: plt.savefig(fname)


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
    if nthreads <= 0: raise ValueError(f'invalid `nthreads`: {nthreads:d}')
    if queue is not None:
      if walltime is None:
        raise ValueError('`walltime` is required when `queue` is set')
      if partition != 'haswell' and partition != 'knl':
        raise ValueError("`partition` must be 'haswell' or 'knl'")
    if (len(seeds) <= 1): raise TypeError('seeds must be a list')
    if None in self._param.values():
      raise ValueError('please set EZmock parameters via `set_param` or `run`')

    jobs = []
    # Generate the script for each seed.
    for s in seeds:
      if s <= 0:
        self._warn(f'omitting non-positive seed: {s:d}')
        continue
      params = copy.deepcopy(self._param)
      params['seed'] = s
      bname = self._get_bname(params)
      odir = f'{self._workdir}/{bname}'
      if not os.path.isdir(odir):
        try: os.mkdir(odir)
        except: raise IOError(f'cannot create directory: {odir}')

      # Check if the job exists already
      par = pkconf = xiconf = bkconf = None
      script = f'{odir}/{self._script}'
      conf = f'{odir}/{self._conf}'
      done = f'{odir}/{self._done}'
      if os.path.isfile(script):
        try:
          with open(conf, 'r') as f:
            par, pkconf, xiconf, bkconf = json.load(f)
          if par == params and os.path.isfile(done) and \
              (clustering == False or (pkconf == self._pkconf and \
              xiconf == self._xiconf and bkconf == self._bkconf)):
            self._warn(f'omitting existing job for seed {s:d}: {script}')
            continue
          if par is not None and par != params:
            os.rename(conf, f'{conf.old}')
            self._warn(('existing EZmock with seed {s:d} will be overwritten, '
                  'due to the change of parameters. The previous settings '
                  f'are moved to {conf}.old'))
        except FileNotFoundError:
          self._warn(('existing EZmock with seed {s:d} detected but the '
                'settings are not found, rerun anyway'))
        if os.path.isfile(done): os.remove(done)

      jobstr = (f'#!/bin/bash\n\nexport OMP_NUM_THREADS={nthreads:d}\n\n'
          f'cd {odir}\n\n')

      run_mock = True
      ofile = f'{odir}/EZmock_{bname}.dat'
      if par == params and os.path.isfile(ofile): run_mock = False
      else: pkconf = xiconf = bkconf = None

      run_rsd = clustering and (self._pkconf['zspace'] or \
          self._xiconf['zspace'] or self._bkconf['zspace'])
      ofile = f'{odir}/EZmock_{bname}_RSD.dat'
      if run_rsd and par == params and os.path.isfile(ofile): run_rsd = False

      if run_mock or run_rsd:
        jobstr += self._mock_cmd(bname, mock=run_mock, rsd=run_rsd,
            params=params)

      if clustering:
        if self._pkconf['rspace']:
          ofile = f'{odir}/PK_EZmock_{bname}.dat'
          if pkconf != self._pkconf or not os.path.isfile(ofile):
            jobstr += self._pk_cmd(rsd=False)
        if self._pkconf['zspace']:
          ofile = f'{odir}/PK_EZmock_{bname}_RSD.dat'
          if pkconf != self._pkconf or not os.path.isfile(ofile):
            jobstr += self._pk_cmd(rsd=True)
        if self._xiconf['rspace']:
          ofile = f'{odir}/2PCF_EZmock_{bname}.dat'
          if xiconf != self._xiconf or not os.path.isfile(ofile):
            jobstr += self._xi_cmd(rsd=False)
        if self._xiconf['zspace']:
          ofile = f'{odir}/2PCF_EZmock_{bname}_RSD.dat'
          if xiconf != self._xiconf or not os.path.isfile(ofile):
            jobstr += self._xi_cmd(rsd=True)
        if self._bkconf['rspace']:
          ofile = f'{odir}/BK_EZmock_{bname}.dat'
          if bkconf != self._bkconf or not os.path.isfile(ofile):
            jobstr += self._bk_cmd(rsd=False)
        if self._bkconf['zspace']:
          ofile = f'{odir}/BK_EZmock_{bname}_RSD.dat'
          if bkconf != self._bkconf or not os.path.isfile(ofile):
            jobstr += self._bk_cmd(rsd=True)

      jobstr += f'echo 1 > {self._done}\n'

      # Save the script and configurations
      with open(script, 'w') as f: f.write(jobstr)
      with open(conf, 'w') as f:
        json.dump([params, self._pkconf, self._xiconf, self._bkconf], f,
            indent=2)
      jobs.append(script)

    # Generate the job script for all realisations.
    njob = len(jobs)
    if njob < 1: raise ValueError('no valid job is found')
    script = f'{self._workdir}/submit_mass_production.sh'
    jobstr = '#!/bin/bash\n'
    if queue is not None:
      jobstr += (f'#SBATCH -n {njob:d}\n#SBATCH -L SCRATCH\n'
        f'#SBATCH -o {self._workdir}/massive_stdout_%j.txt\n'
        f'#SBATCH -e {self._workdir}/massive_stderr_%j.txt\n')
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


  def restore(self):
    """
    Restore the parameters of previous runs from existing files.
    """
    from glob import glob
    paths = glob(f'{self._workdir}/B*G*Z*N*_b*d*r*c*_seed*')
    for p in paths:
      if os.path.isfile(f'{p}/{self._done}'):
        with open(f'{p}/{self._conf}', 'r') as f:
          par, _, _, _ = json.load(f)
        if not par in self._history: self._history.append(par)


  def params(self):
    """
    Print the current set of EZmock parameters.
    """
    for key, value in self._param.items(): print(key, '=', value)


  def history(self):
    """
    Print the histories of the EZmock parameter sets.
    """
    for i, param in enumerate(self._history): print(f'{i:d}:', param)


  def clear(self, slicer):
    """
    Clear history entries defined by `slicer`.

    Parameters
    ----------
    slicer:
        The slice of the histories to be cleared.
        It must be generated using the `slice` function.
    """
    if type(slicer).__name__ != 'slice':
      raise TypeError('slicer must be generated using the `slice` function')
    del(self._history[slicer])


  @property
  def workdir(self):
    """
    Working directory for generating EZmocks.
    """
    return self._workdir

  @workdir.setter
  def workdir(self, path):
    if path is None:
      raise ValueError('working directory not set')
    if not os.path.isdir(path):
      try: os.mkdir(path)
      except: raise IOError(f'cannot create working directory: {path}')
    self._workdir = os.path.abspath(path)


  @property
  def exe(self):
    """
    Path of the EZmock executable.
    """
    return self._ez_exe

  @exe.setter
  def exe(self, path):
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
      raise IOError(f'invalid EZmock executable: {path}')
    self._ez_exe = os.path.abspath(path)


  @property
  def pk_exe(self):
    """
    Path of the powspec executable.
    """
    return self._pk_exe

  @pk_exe.setter
  def pk_exe(self, path):
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
      self._warn(f'Power spectrum disabled due to invalid executable: {path}')
      self._pk_exe = None
    else: self._pk_exe = os.path.abspath(path)


  @property
  def xi_exe(self):
    """
    Path of the FCFC executable for periodic boxes.
    """
    return self._xi_exe

  @xi_exe.setter
  def xi_exe(self, path):
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
      self._warn(f'2PCF disabled due to invalid executable: {path}')
      self._xi_exe = None
    else: self._xi_exe = os.path.abspath(path)


  @property
  def bk_exe(self):
    """
    Path of the bispec executable.
    """
    return self._bk_exe

  @bk_exe.setter
  def bk_exe(self, path):
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
      self._warn(f'Bispectrum disabled due to invalid executable: {path}')
      self._bk_exe = None
    else: self._bk_exe = os.path.abspath(path)


  def _get_bname(self, param):
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


  def _mock_cmd(self, bname, mock=True, rsd=True, params=None):
    """
    Generate the command for constructing the EZmock catalogue.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.
    mock: bool optional
        Indicate whether to run EZmock.
    rsd: bool, optional
        Indicate whether to apply redshift space distortions.
    params: dict, optional
        A dictionary of the EZmock parameters.

    Returns
    -------
    The command as a string
    """
    from cosmoEZ import flatLCDM

    if params is None: params = self._param
    jobstr = ''

    # Compute structure growth parameters
    cosmo = flatLCDM(omega_m = params['omega_m'])
    z1 = 1 + params['redshift']
    a = 1. / z1
    a_init = 1. / (1 + params['z_init'])
    grow2z0 = cosmo.growth2(a, a_init=a_init)
    hubble = cosmo.hubble(a)
    zdist = cosmo.growthf(a) * hubble * a

    # Generate the command for running EZmock
    if mock:
      jobstr += ("echo '&EZmock_v0_input\n"
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
          f" | {self._ez_exe} || exit\n\n")

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


  def _pk_cmd(self, bname=None, rsd=True):
    """
    Generate the command for computing power spectrum of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.
    rsd: bool, optional
        Indicate whether to compute the redshift-space power spectrum.

    Returns
    -------
    The command as a string
    """
    if bname is None: bname = self._bname

    if rsd:
      poles = '[' + ','.join(f'{i:d}' for i in self._pkconf['zell']) + ']'
      ifile = f'EZmock_{bname}_RSD.dat'
    else:
      poles = '[' + ','.join(f'{i:d}' for i in self._pkconf['rell']) + ']'
      ifile = f'EZmock_{bname}.dat'
    ofile = f'PK_{ifile}'
    bsize = self._param['boxsize']

    jobstr = (f'{self._pk_exe} -d {ifile} --data-formatter "%lf %lf %lf" '
        f"-p '[($1+{bsize:g})%{bsize:g},($2+{bsize:g})%{bsize:g},"
        f"($3+{bsize:g})%{bsize:g}]' -s T -B {bsize:g} "
        f"-G {self._pkconf['ngrid']:d} -n 1 -i F -l '{poles}' -k 0 "
        f"-K {self._pkconf['kmax']:g} -b {self._pkconf['dk']:g} "
        f"-a {ofile} || exit\n\n")

    return jobstr


  def _xi_cmd(self, bname=None, rsd=True):
    """
    Generate the command for computing 2-point correlation function of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.
    rsd: bool, optional
        Indicate whether to compute the redshift-space power spectrum.

    Returns
    -------
    The command as a string
    """
    if bname is None: bname = self._bname
    if rsd:
      poles = '[' + ','.join(f'{i:d}' for i in self._xiconf['zell']) + ']'
      ifile = f'EZmock_{bname}_RSD.dat'
    else:
      poles = '[' + ','.join(f'{i:d}' for i in self._xiconf['rell']) + ']'
      ifile = f'EZmock_{bname}.dat'
    ofile = f'2PCF_{ifile}'
    bsize = self._param['boxsize']

    jobstr = (f"{self._xi_exe} -i {ifile} -l D -f '%lf %lf %lf' "
        f"-x '[$1,$2,$3]' -b {bsize:g} -B 1 -p DD -P {ofile}.dd -e 'DD/@@-1' "
        f"-E {ofile}.xi2d -m '{poles}' -M {ofile} --s-min 0 "
        f"--s-max {self._xiconf['rmax']:g} --s-step {self._xiconf['dr']:g} "
        f"--mu-num {self._xiconf['nmu']:d} --dist-prec 0 -S 0 || exit\n\n")

    return jobstr


  def _bk_cmd(self, bname=None, rsd=True):
    """
    Generate the command for computing power spectrum of EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.
    rsd: bool, optional
        Indicate whether to compute the redshift-space power spectrum.

    Returns
    -------
    The command as a string
    """
    if bname is None: bname = self._bname
    if rsd: ifile = f'EZmock_{bname}_RSD.dat'
    else: ifile = f'EZmock_{bname}.dat'
    ofile = f'BK_{ifile}'
    bsize = self._param['boxsize']

    jobstr = (f'{self._bk_exe} -i {ifile} -s 7 --x-min 0 --y-min 0 --z-min 0 '
        f'--x-max {bsize:g} --y-max {bsize:g} --z-max {bsize:g} -b 0'
        f" -B {bsize:g} -g {self._bkconf['ngrid']:d} -w 1 -x 0 "
        f"-p {self._bkconf['k1'][0]:g} -P {self._bkconf['k1'][1]:g} "
        f"-q {self._bkconf['k2'][0]:g} -Q {self._bkconf['k2'][1]:g} "
        f"-n {self._bkconf['nbin']:d} -o {ofile} -y || exit\n\n")

    return jobstr


  def _warn(self, message):
    """
    Print the warning message.

    Parameters
    ----------
    message: str
        The message to be printed.
    """
    print('\x1b[31;1mWarning:\x1b[0m ' + message)

