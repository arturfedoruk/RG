# ``betas.m``
Contains beta-functions, up to |4,4,4|3,3|3| loop configuration and declarations of ``_dot[loops_]`` functions.

# ``phase_spases.nb``
Notebook with picture of different phase spaces

# ``loops_compare.nb``
Notebook where contributions of loops in different variables are compared

# ``fit_tables.nb``
Notebook for effortless (almost) transforming fit results ``../fit_outputs/``, obtained by ``../fit.cpp``, into TeX for for the paper. Just insert 
"For table" entries of the form ``{x, y, z, ...}`` from ``output_tree___.txt`` into the ``_Fit`` wolfram tables and run the cell.

# ``parametric_N_theoretical.nb``
Notebook, where parametric and theoretical uncertainties of parameters of lambda running (interception with lambda=0, scale mu:beta=0 and value lambda(mu:beta=0))
are calculated. Parametric uncertainties are calculated by generating random initial values at mu=173, according to fits and propogation of errors computed
by ``../fit.cpp``, and then calculating the listed values for the running curve. Values are then stored at ``sols_data_sofar.dat`` in order not to generate it
multiple times. To get started, insert central values of MSbar into ``means[loops_]`` and their uncertainties into ``sigmas[loops_]``. 
The former is obtained from the first number in ``For table`` lists in ``output_tree___.txt``, the latter - from propogation of errors ``Delta ai`` from the same file.
Theoretical uncertainty is calculated by running values, obtained by ``../theor_err.cpp``, calculating parameters of lambda running, and
picking maximal and minimal values. The notebook is pretty large as a result of a need to graph both theoretical and parametric at the same plot (defenately could be done more elegantly)
