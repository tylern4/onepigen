onepigen is a single charged pion event generator which supports cross
sections generated via AO, daresbury, or MAID tables, in addition to
supporting Fermi motion in the target as well as generation of a
spectator particle.  In addition, it might be modified in the future
to support arbitrary cross sections via shaping an interpolated model
of a sampled cross section.

The model used for Fermi motion supports many settings.  First, the
energy component can be determined in a variety of ways, even though
the off-shellness in the energy component is not significant in that
many different reasonable estimates yield indistinguishable results,
at least to the level of precision currently available in cross
section measurements.  Second, the Hulthen, Paris, and CD-Bonn
momentum distributions can be used to model the Fermi momentum
distribution.  It is not clear to me whether or not the asymmetrical
angular states should be handled or if they have any meaningful
effects.

The program is controlled via the standard input, and accepts
Fortran's namelist I/O formatting.  This seems to be the most
efficient way to handle the large number of control parameters the
program accepts as opposed to flags supplied at the command line.

## Build and run

Clone the repo and build onepigen
```
git clone https://github.com/tylern4/onepigen.git
cd onepigen
make
```

Inside the onepigen directory add the executable `onepigen_lund` to your path so the python script can find it.
```
export PATH=${PATH}:${PWD}/bin
```

Run the python script
```
./onepigen --trig=100 --experiment=rgb
```


All options
```
usage: onepigen [-h] [--trig TRIG] [--docker] [--experiment EXPERIMENT] [--flag_ehel FLAG_EHEL] [--epirea EPIREA] [--radeffects RADEFFECTS] [--mm_cut MM_CUT] [--reg1 REG1] [--reg2 REG2]
                [--reg3 REG3] [--reg4 REG4] [--t_targ T_TARG] [--r_targ R_TARG] [--vertex_x VERTEX_X] [--vertex_y VERTEX_Y] [--vz VZ] [--beam_energy BEAM_ENERGY] [--q2_min Q2_MIN]
                [--q2_max Q2_MAX] [--w_min W_MIN] [--w_max W_MAX] [--ep_min EP_MIN] [--ep_max EP_MAX] [--delta DELTA] [--fmcall FMCALL] [--sigr_max SIGR_MAX] [--fermi_motion FERMI_MOTION]
                [--fermi_energy FERMI_ENERGY]

Run onepigen

options:
  -h, --help            show this help message and exit
  --trig TRIG
  --docker
  --experiment EXPERIMENT
  --flag_ehel FLAG_EHEL
  --epirea EPIREA
  --radeffects RADEFFECTS
  --mm_cut MM_CUT
  --reg1 REG1
  --reg2 REG2
  --reg3 REG3
  --reg4 REG4
  --t_targ T_TARG
  --r_targ R_TARG
  --vertex_x VERTEX_X
  --vertex_y VERTEX_Y
  --vz VZ
  --beam_energy BEAM_ENERGY
  --q2_min Q2_MIN
  --q2_max Q2_MAX
  --w_min W_MIN
  --w_max W_MAX
  --ep_min EP_MIN
  --ep_max EP_MAX
  --delta DELTA
  --fmcall FMCALL
  --sigr_max SIGR_MAX
  --fermi_motion FERMI_MOTION
  --fermi_energy FERMI_ENERGY
```
