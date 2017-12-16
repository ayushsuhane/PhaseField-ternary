import os
delete = 'rm outphase_*.txt line_*.txt phase_*.txt comp_*.txt grain_*.txt'
os.system(delete)
#c = 'gcc phasefield.c readinputs.c initialize.c allocate.c eqmole.c mucalc.c convert.c uniformrn.c setup.c visual.c solver.c superimpose.c -lgsl -lgslcblas -lm'
#c = 'gcc -Wall -g -pg phasefield.c readinputs.c initialize.c constant_eval.c eval_free_energy.c pd_solve.c eqmole.c mucalc.c convert.c allocate.c uniformrn.c setup.c visual.c superimpose.c solver.c -lgsl -lgslcblas -lm -O2 -O3'
#c = 'gcc -Wall -g -pg phasefield.c readinputs.c initialize.c initialize_fn.c solver_functions.c allocate.c uniformrn.c stencil.c setup.c visual.c superimpose.c solver.c -lgsl -lgslcblas -lm -O2 -O3'
#c = 'gcc -Wall -g -w -pg phasefield.c readinputs.c initialize.c initialize_fn.c allocate.c solver_functions.c uniformrn.c  setup.c visual.c stencil.c superimpose.c solver.c -lgsl -lgslcblas -lm -O2 -O3'
#c = 'gcc -Wall -g -w -pg phasefield.c readinputs.c initialize.c initialize_fn.c allocate.c setup.c uniformrn.c solver.c superimpose.c visual.c stencil.c solver_functions.c check_phases.c -lgsl -lgslcblas -lm -O2 -O3'
#c = 'gcc -Wall -g -w -pg phasefield.c readinputs.c allocate.c setup_fn.c setup.c uniformrn.c stencil.c superimpose.c visual.c solver.c solver_functions.c memory.c -lgsl -lgslcblas -lm -O2 -O3 -funroll-loops -o PHASEFIELD'
c = 'gcc -Wall -w  phasefield.c readinputs.c allocate.c setup_fn.c setup.c uniformrn.c stencil.c solver_setup.c visual.c solver.c solver_functions.c memory.c -lgsl -lgslcblas -lm -O2 -O3 -funroll-loops -o PHASEFIELD'

os.system(c)
