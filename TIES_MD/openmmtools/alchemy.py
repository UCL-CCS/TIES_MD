import copy
import itertools
import logging

import openmmtools
from openmmtools.constants import ONE_4PI_EPS0
from openmmtools import states
try:
    import openmm
    from openmm import unit
except ImportError:  # OpenMM < 7.6
    from simtk import openmm, unit


logger = logging.getLogger(__name__)

class AlchemicalStateError(states.GlobalParameterError):
    """Error raised by an AlchemicalState."""
    pass

class ModifiedAlchemicalState(states.GlobalParameterState):
    """Represent an alchemical state.
    The alchemical parameters modify the Hamiltonian and affect the
    computation of the energy. Alchemical parameters that have value
    None are considered undefined, which means that applying this
    state to System and Context that have that parameter as a global
    variable will raise an AlchemicalStateError.
    Parameters
    ----------
    parameters_name_suffix : str, optional
        If specified, the state will control a modified version of the global
        parameters with the name ``parameter_name + '_' + parameters_name_suffix``.
        When this is the case, the normal parameters are not accessible.
    lambda_sterics : float, optional
        Scaling factor for ligand sterics (Lennard-Jones and Halgren)
        interactions (default is 1.0).
    lambda_electrostatics : float, optional
        Scaling factor for ligand charges, intrinsic Born radii, and surface
        area term (default is 1.0).
    lambda_bonds : float, optional
        Scaling factor for alchemically-softened bonds (default is 1.0).
    lambda_angles : float, optional
        Scaling factor for alchemically-softened angles (default is 1.0).
    lambda_torsions : float, optional
        Scaling factor for alchemically-softened torsions (default is 1.0).
    Attributes
    ----------
    lambda_sterics
    lambda_electrostatics
    lambda_bonds
    lambda_angles
    lambda_torsions
    Examples
    --------
    Create an alchemically modified system.
    >>> from openmmtools import testsystems
    >>> factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
    >>> alanine_vacuum = testsystems.AlanineDipeptideVacuum().system
    >>> alchemical_region = AlchemicalRegion(alchemical_atoms=range(22))
    >>> alanine_alchemical_system = factory.create_alchemical_system(reference_system=alanine_vacuum,
    ...                                                              alchemical_regions=alchemical_region)
    Create a completely undefined alchemical state.
    >>> alchemical_state = AlchemicalState()
    >>> print(alchemical_state.lambda_sterics)
    None
    >>> alchemical_state.apply_to_system(alanine_alchemical_system)
    Traceback (most recent call last):
    ...
    openmmtools.alchemy.AlchemicalStateError: The system parameter lambda_electrostatics is not defined in this state.
    Create an AlchemicalState that matches the parameters defined in
    the System.
    >>> alchemical_state = AlchemicalState.from_system(alanine_alchemical_system)
    >>> alchemical_state.lambda_sterics
    1.0
    >>> alchemical_state.lambda_electrostatics
    1.0
    >>> print(alchemical_state.lambda_angles)
    None
    AlchemicalState implement the IComposableState interface, so it can be
    used with CompoundThermodynamicState. All the alchemical parameters are
    accessible through the compound state.
    >>> import openmm
    >>> from openmm import unit
    >>> thermodynamic_state = states.ThermodynamicState(system=alanine_alchemical_system,
    ...                                                 temperature=300*unit.kelvin)
    >>> compound_state = states.CompoundThermodynamicState(thermodynamic_state=thermodynamic_state,
    ...                                                    composable_states=[alchemical_state])
    >>> compound_state.lambda_sterics
    1.0
    You can control the parameters in the OpenMM Context in this state by
    setting the state attributes.
    >>> compound_state.lambda_sterics = 0.5
    >>> integrator = openmm.VerletIntegrator(1.0*unit.femtosecond)
    >>> context = compound_state.create_context(integrator)
    >>> context.getParameter('lambda_sterics')
    0.5
    >>> compound_state.lambda_sterics = 1.0
    >>> compound_state.apply_to_context(context)
    >>> context.getParameter('lambda_sterics')
    1.0
    You can express the alchemical parameters as a mathematical expression
    involving alchemical variables. Here is an example for a two-stage function.
    >>> compound_state.set_alchemical_variable('lambda', 1.0)
    >>> compound_state.lambda_sterics = AlchemicalFunction('step_hm(lambda - 0.5) + 2*lambda * step_hm(0.5 - lambda)')
    >>> compound_state.lambda_electrostatics = AlchemicalFunction('2*(lambda - 0.5) * step(lambda - 0.5)')
    >>> for l in [0.0, 0.25, 0.5, 0.75, 1.0]:
    ...     compound_state.set_alchemical_variable('lambda', l)
    ...     print(compound_state.lambda_sterics)
    0.0
    0.5
    1.0
    1.0
    1.0
    """

    _GLOBAL_PARAMETER_ERROR = AlchemicalStateError

    # -------------------------------------------------------------------------
    # Lambda properties
    # -------------------------------------------------------------------------

    class _LambdaParameter(states.GlobalParameterState.GlobalParameter):
        """A global parameter in the interval [0, 1] with standard value 1."""

        def __init__(self, parameter_name):
            super().__init__(parameter_name, standard_value=1.0, validator=self.lambda_validator)

        @staticmethod
        def lambda_validator(self, instance, parameter_value):
            if parameter_value is None:
                return parameter_value
            if not (0.0 <= parameter_value <= 1.0):
                print('{} must be between 0 and 1.'.format(self.parameter_name))
            return float(parameter_value)

    lambda_sterics = _LambdaParameter('lambda_sterics')
    lambda_electrostatics = _LambdaParameter('lambda_electrostatics')
    lambda_bonds = _LambdaParameter('lambda_bonds')
    lambda_angles = _LambdaParameter('lambda_angles')
    lambda_torsions = _LambdaParameter('lambda_torsions')

    @classmethod
    def from_system(cls, system, *args, **kwargs):
        """Constructor reading the state from an alchemical system.
        Parameters
        ----------
        system : openmm.System
            An alchemically modified system in a defined alchemical state.
        parameters_name_suffix : str, optional
            If specified, the state will search for a modified
            version of the alchemical parameters with the name
            ``parameter_name + '_' + parameters_name_suffix``.
        Returns
        -------
        The AlchemicalState object representing the alchemical state of
        the system.
        Raises
        ------
        AlchemicalStateError
            If the same parameter has different values in the system, or
            if the system has no lambda parameters.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        return super().from_system(system, *args, **kwargs)

    def set_alchemical_parameters(self, new_value):
        """Set all defined lambda parameters to the given value.
        The undefined parameters (i.e. those being set to None) remain
        undefined.
        Parameters
        ----------
        new_value : float
            The new value for all defined parameters.
        """
        for parameter_name in self._parameters:
            if self._parameters[parameter_name] is not None:
                setattr(self, parameter_name, new_value)

    # -------------------------------------------------------------------------
    # Function variables
    # -------------------------------------------------------------------------

    def get_function_variable(self, variable_name):
        """Return the value of the function variable.
        Function variables are variables entering mathematical expressions
        specified with ``AlchemicalFunction``, which can be use to enslave
        a lambda parameter to arbitrary variables.
        Parameters
        ----------
        variable_name : str
            The name of the function variable.
        Returns
        -------
        variable_value : float
            The value of the function variable.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        return super().get_function_variable(variable_name)

    def set_function_variable(self, variable_name, new_value):
        """Set the value of the function variable.
        Function variables are variables entering mathematical expressions
        specified with ``AlchemicalFunction``, which can be use to enslave
        a lambda parameter to arbitrary variables.
        Parameters
        ----------
        variable_name : str
            The name of the function variable.
        new_value : float
            The new value for the variable.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        super().set_function_variable(variable_name, new_value)

    def get_alchemical_variable(self, variable_name):
        """Return the value of the alchemical parameter.
        .. warning:
            This is deprecated. Use ``get_function_variable`` instead.
        Parameters
        ----------
        variable_name : str
            The name of the alchemical variable.
        Returns
        -------
        variable_value : float
            The value of the alchemical variable.
        """
        import warnings
        warnings.warn('AlchemicalState.get_alchemical_variable is deprecated. '
                      'Use AlchemicalState.get_function_variable instead.')
        return super().get_function_variable(variable_name)

    def set_alchemical_variable(self, variable_name, new_value):
        """Set the value of the alchemical variable.
        .. warning:
            This is deprecated. Use ``set_function_variable`` instead.
        Parameters
        ----------
        variable_name : str
            The name of the alchemical variable.
        new_value : float
            The new value for the variable.
        """
        import warnings
        warnings.warn('AlchemicalState.get_alchemical_variable is deprecated. '
                      'Use AlchemicalState.get_function_variable instead.')
        super().set_function_variable(variable_name, new_value)

    # -------------------------------------------------------------------------
    # IComposableState interface
    # -------------------------------------------------------------------------

    def apply_to_system(self, system):
        """Set the alchemical state of the system to this.
        Parameters
        ----------
        system : openmm.System
            The system to modify.
        Raises
        ------
        AlchemicalStateError
            If the system does not have the required lambda global variables.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        super().apply_to_system(system)

    def check_system_consistency(self, system):
        """Check if the system is in this alchemical state.
        It raises a AlchemicalStateError if the system is not consistent
        with the alchemical state.
        Parameters
        ----------
        system : openmm.System
            The system to test.
        Raises
        ------
        AlchemicalStateError
            If the system is not consistent with this state.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        super().check_system_consistency(system)

    def apply_to_context(self, context):
        """Put the Context into this AlchemicalState.
        Parameters
        ----------
        context : openmm.Context
            The context to set.
        Raises
        ------
        AlchemicalStateError
            If the context does not have the required lambda global variables.
        """
        # The function is redefined here only to provide more specific documentation for this method.
        super().apply_to_context(context)

class ModifiedAbsoluteAlchemicalFactory(openmmtools.alchemy.AbsoluteAlchemicalFactory):
    """super?
        Please see ...
    """

    def _get_sterics_energy_expressions(self, lambda_variable_suffixes, softcore):
        """Return the energy expressions for sterics.

        Parameters
        ----------
        lambda_variable_suffixes : List[str]
            A list with suffixes for the global variable "lambda_sterics" that
            will control the energy. If no suffix is necessary (i.e. there are
            no multiple alchemical regions) just set lambda_variable_suffixes[0] = ''.
            If the list has more than one element, the energy is controlled by
            the multiplication of lambda_sterics_suffix1 * lambda_sterics_suffix2.
        softcore : dict
            A dictionary which contains softcore parameters as floats
        """

        # Sterics mixing rules.
        if lambda_variable_suffixes[0] == '':
            lambda_variable_name = 'lambda_sterics'
        else:
            if len(lambda_variable_suffixes) > 1:
                lambda_variable_name = 'lambda_sterics{0}*lambda_sterics{1}'.format(
                    lambda_variable_suffixes[0], lambda_variable_suffixes[1])
            else:
                lambda_variable_name = 'lambda_sterics{}'.format(lambda_variable_suffixes[0])

        sterics_mixing_rules = ('epsilon = sqrt(epsilon1*epsilon2);'  # Mixing rule for epsilon.
                                'sigma = 0.5*(sigma1 + sigma2);')  # Mixing rule for sigma.

        # Soft-core Lennard-Jones.
        exceptions_sterics_energy_expression = ('U_sterics;'
                                                'U_sterics = (({0})^{1})*4*epsilon*x*(x-1.0);'
                                                'x = (sigma/reff_sterics)^6;'
                                                # Effective softcore distance for sterics.
                                                'reff_sterics = sigma*(({2}*(1.0-({0}))^{3} + (r/sigma)^{4}))^(1/{4});')\
                                                .format(lambda_variable_name, softcore['a'], softcore['alpha'], softcore['b'], softcore['c'])
        # Define energy expression for electrostatics.
        return sterics_mixing_rules, exceptions_sterics_energy_expression

    def _get_electrostatics_energy_expressions(self, reference_force, lambda_variable_suffixes, softcore):
        """Return the energy expressions for electrostatics.

        This private function assumes self._alchemical_pme_treatment != 'exact'
        as there's no electrostatics CustomNondondedForce in this case, and
        lambda_electrostatics is modeled through an offset parameter in a
        NonbondedForce.

        Parameters
        ----------
        lambda_variable_suffixes : List[str]
            A list with suffixes for the global variable "lambda_electrostatics" that
            will control the energy. If no suffix is necessary (i.e. there are
            no multiple alchemical regions) just set lambda_variable_suffixes[0] = ''.
            If the list has more than one element, the energy is controlled by
            the multiplication of lambda_electrostatics_suffix1 * lambda_electrostatics_suffix2.
        softcore : dict
            A dictionary which contains softcore parameters as floats
        """
        if lambda_variable_suffixes[0] == '':
            lambda_variable_name = 'lambda_electrostatics'
        else:
            if len(lambda_variable_suffixes) > 1:
                lambda_variable_name = 'lambda_electrostatics{0}*lambda_electrostatics{1}'.format(
                    lambda_variable_suffixes[0], lambda_variable_suffixes[1])
            else:
                lambda_variable_name = 'lambda_electrostatics{}'.format(lambda_variable_suffixes[0])

        # The final expression will be prefix + method + suffix.
        electrostatics_prefix = ('U_electrostatics;U_electrostatics=(({0})^{1})*ONE_4PI_EPS0*chargeprod')\
            .format(lambda_variable_name, softcore['d'])

        # Effective softcore distance for electrostatics (common to all methods).
        electrostatics_suffix = ('reff_electrostatics = sigma*(({0}*(1.0-({1}))^{2} + (r/sigma)^{3}))^(1/{3});'
                                 ' ONE_4PI_EPS0 = {4};').format(softcore['beta'], lambda_variable_name, softcore['e'],
                                                                softcore['f'], ONE_4PI_EPS0)  # Already in OpenMM units.

        # Define mixing rules.
        electrostatics_mixing_rules = ('chargeprod = charge1*charge2;'  # Mixing rule for charges.
                                       'sigma = 0.5*(sigma1 + sigma2);')  # Mixing rule for sigma.

        # Standard Coulomb expression with softened core. This is used
        #   - When the nonbonded method of the reference force is NoCutoff.
        #   - When alchemical_pme_treatment is set to 'coulomb'.
        #   - With 1-4 exceptions, unless self.consistent_exceptions is True.
        coulomb_expression = '/reff_electrostatics;'

        # Select electrostatics functional form based on nonbonded method.
        nonbonded_method = reference_force.getNonbondedMethod()

        # Soft-core Coulomb.
        if nonbonded_method in [openmm.NonbondedForce.NoCutoff]:
            electrostatics_method_expression = coulomb_expression
        # Reaction-field electrostatics.
        elif nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic]:
            electrostatics_method_expression = self._get_reaction_field_unique_expression(reference_force)
        # PME electrostatics.
        elif nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
            # Ewald direct-space electrostatics.
            if self.alchemical_pme_treatment == 'direct-space':
                electrostatics_method_expression = self._get_pme_direct_space_unique_expression(reference_force)
            # Use switched standard Coulomb potential, following MTS scheme described in
            # http://dx.doi.org/10.1063/1.1385159
            elif self.alchemical_pme_treatment == 'coulomb':
                electrostatics_method_expression = coulomb_expression
            else:
                raise ValueError("Unknown alchemical_pme_treatment scheme '{}'".format(self.alchemical_pme_treatment))
        else:
            raise ValueError("Nonbonded method {} not supported yet.".format(nonbonded_method))

        # Define energy expression for 1,4 electrostatic exceptions.
        exceptions_electrostatics_energy_expression = electrostatics_prefix
        if self.consistent_exceptions:
            exceptions_electrostatics_energy_expression += electrostatics_method_expression
        else:
            exceptions_electrostatics_energy_expression += coulomb_expression
        exceptions_electrostatics_energy_expression += electrostatics_suffix

        # Define energy expression for electrostatics.
        electrostatics_energy_expression = (electrostatics_prefix + electrostatics_method_expression +
                                            electrostatics_suffix + electrostatics_mixing_rules)

        return electrostatics_energy_expression, exceptions_electrostatics_energy_expression

    def _alchemically_modify_NonbondedForce(self, reference_force, alchemical_regions, alchemical_regions_interactions):
        """Create alchemically-modified version of NonbondedForce.

        Parameters
        ----------
        reference_force : openmm.NonbondedForce
            The reference NonbondedForce to be alchemically modify.
        alchemical_region : AlchemicalRegion
            The alchemical region containing the indices of the atoms to
            alchemically modify.
        alchemical_regions_interactions : Set[Tuple[int, int]], optional
            Set of alchemical region index pairs for interacting regions.
            By default, all alchemical regions interact only with the
            non-alchemical environment.

        Returns
        -------
        nonbonded_force : openmm.NonbondedForce
            The force responsible for interactions and exceptions of non-alchemical atoms.
        aa_sterics_custom_nonbonded_force : openmm.CustomNonbondedForce
            The force responsible for sterics interactions of alchemical/alchemical atoms.
        aa_electrostatics_custom_nonbonded_force : openmm.CustomNonbondedForce
            The force responsible for electrostatics interactions of alchemical/alchemical
            atoms.
        na_sterics_custom_nonbonded_force : openmm.CustomNonbondedForce
            The force responsible for sterics interactions of non-alchemical/alchemical atoms.
        na_electrostatics_custom_nonbonded_force : openmm.CustomNonbondedForce
            The force responsible for electrostatics interactions of non-alchemical/alchemical
            atoms.
        aa_sterics_custom_bond_force : openmm.CustomBondForce
            The force responsible for sterics exceptions of alchemical/alchemical atoms.
        aa_electrostatics_custom_bond_force : openmm.CustomBondForce
            The force responsible for electrostatics exceptions of alchemical/alchemical
            atoms.
        na_sterics_custom_bond_force : openmm.CustomBondForce
            The force responsible for sterics exceptions of non-alchemical/alchemical atoms.
        na_electrostatics_custom_bond_force : openmm.CustomBondForce
            The force responsible for electrostatics exceptions of non-alchemical/alchemical
            atoms.

        References
        ----------
        [1] Pham TT and Shirts MR. Identifying low variance pathways for free
        energy calculations of molecular transformations in solution phase.
        JCP 135:034114, 2011. http://dx.doi.org/10.1063/1.3607597

        """
        # TODO Change softcore_beta to a dimensionless scalar to multiply some intrinsic length-scale, like Lennard-Jones alpha.
        # TODO Try using a single, common "reff" effective softcore distance for both Lennard-Jones and Coulomb.

        forces_by_lambda = {}
        all_alchemical_atoms = set()

        # Don't create a force if there are no alchemical atoms.
        for alchemical_region in alchemical_regions:
            if not len(alchemical_region.alchemical_atoms) == 0:
                all_alchemical_atoms.update(alchemical_region.alchemical_atoms)

        # Don't create a force if there are no alchemical atoms.
        if len(all_alchemical_atoms) == 0:
            return {'': [copy.deepcopy(reference_force)]}

        # Create a set of all the non-alchemical atoms only.
        all_atomset = set(range(reference_force.getNumParticles()))
        nonalchemical_atomset = all_atomset.difference(all_alchemical_atoms)

        # -------------------------------------------------------------
        # Perform tasks that do not need to be repeated for all regions.
        # -------------------------------------------------------------

        # Define energy expression for electrostatics based on nonbonded method.
        nonbonded_method = reference_force.getNonbondedMethod()
        is_ewald_method = nonbonded_method in [openmm.NonbondedForce.Ewald,
                                               openmm.NonbondedForce.PME]
        is_rf_method = nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic,
                                            openmm.NonbondedForce.CutoffNonPeriodic]
        is_periodic_method = is_ewald_method or nonbonded_method == openmm.NonbondedForce.CutoffPeriodic
        use_exact_pme_treatment = is_ewald_method and self.alchemical_pme_treatment == 'exact'

        # Warn about reaction field.
        if is_rf_method:
            logger.warning('Reaction field support is still experimental. For free energy '
                           'calculations in explicit solvent, we suggest using PME for now.')

        # Check that PME treatment is supported with the region's parameters.
        if use_exact_pme_treatment:
            for alchemical_region in alchemical_regions:
                err_msg = ' not supported with exact treatment of Ewald electrostatics.'
                if not alchemical_region.annihilate_electrostatics:
                    raise ValueError('Decoupled electrostatics is' + err_msg)
                if self.consistent_exceptions:
                    raise ValueError('Consistent exceptions are' + err_msg)
                if (alchemical_region.softcore_beta, alchemical_region.softcore_d, alchemical_region.softcore_e) != (0, 1, 1):
                    raise ValueError('Softcore electrostatics is' + err_msg)

        # Create a copy of the NonbondedForce to handle particle interactions and
        # 1,4 exceptions between non-alchemical/non-alchemical atoms (nn).
        nonbonded_force = copy.deepcopy(reference_force)

        # Fix any issues in reference force with Lennard-Jones sigma = 0 (epsilon = 0),
        # which should have sigma > 0.
        for particle_index in range(reference_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = reference_force.getParticleParameters(particle_index)
            # Check particle sigma is not zero.
            if sigma == 0.0 * unit.angstrom:
                warning_msg = 'particle %d has Lennard-Jones sigma = 0 (charge=%s, sigma=%s, epsilon=%s); setting sigma=1A'
                logger.warning(warning_msg % (particle_index, str(charge), str(sigma), str(epsilon)))
                sigma = 1.0 * unit.angstrom
                # Fix it.
                nonbonded_force.setParticleParameters(particle_index, charge, sigma, epsilon)

        # Same for the exceptions.
        for exception_index in range(reference_force.getNumExceptions()):
            # Retrieve parameters.
            [iatom, jatom, chargeprod, sigma, epsilon] = reference_force.getExceptionParameters(exception_index)
            # Check particle sigma is not zero.
            if sigma == 0.0 * unit.angstrom:
                warning_msg = 'exception %d has Lennard-Jones sigma = 0 (iatom=%d, jatom=%d, chargeprod=%s, sigma=%s, epsilon=%s); setting sigma=1A'
                logger.warning(warning_msg % (exception_index, iatom, jatom, str(chargeprod), str(sigma), str(epsilon)))
                sigma = 1.0 * unit.angstrom
                # Fix it.
                nonbonded_force.setExceptionParameters(exception_index, iatom, jatom, chargeprod, sigma, epsilon)

        if use_exact_pme_treatment:
            # Exclude noninteracting alchemical regions from seeing each other in the nonbonded
            for x, y in itertools.combinations(range(len(alchemical_regions)), 2):
                if (x, y) not in alchemical_regions_interactions:
                    for atom1 in alchemical_regions[x].alchemical_atoms:
                        for atom2 in alchemical_regions[y].alchemical_atoms:
                            nonbonded_force.addException(atom1, atom2, 0.0, 1.0, 0.0, True)
                else:
                    region_names = (alchemical_regions[x].name, alchemical_regions[y].name)
                    logger.debug(f'Adding a exact PME electrostatic interaction group between groups {region_names}.')
                    del region_names

            # With exact PME treatment, particle electrostatics is handled through offset parameters.
            for alchemical_region in alchemical_regions:
                if alchemical_region.name is None:
                    nonbonded_force.addGlobalParameter('lambda_electrostatics', 1.0)
                else:
                    nonbonded_force.addGlobalParameter(f'lambda_electrostatics_{alchemical_region.name}', 1.0)

        # Make of list of all single and double permutations of alchemical regions.
        single_regions = [[alchemical_region] for alchemical_region in alchemical_regions]
        if len(alchemical_regions_interactions) == 0:
            pair_regions = []
        else:
            # Only generate pairs of alchemical regions specified by alchemical_regions_interactions.
            pair_regions = [[alchemical_regions[x[0]], alchemical_regions[x[1]]] for x in
                            alchemical_regions_interactions]

        # Iterate over all single and double permutations of alchemical regions to build all interactions.
        for alchemical_regions_pairs in single_regions+pair_regions:

            # Make a list of region names for the alchemical regions interactions which are being built.
            lambda_var_suffixes = []
            for alchemical_region in alchemical_regions_pairs:
                if alchemical_region.name is None:
                    lambda_var_suffixes.append('')
                else:
                    lambda_var_suffixes.append('_' + alchemical_region.name)

            # --------------------------------------------------
            # Determine energy expression for all custom forces
            # --------------------------------------------------

            #Assumes softcore params the same between regions
            softcore_dict = {'alpha': alchemical_regions_pairs[0].softcore_alpha, 'beta': alchemical_regions_pairs[0].softcore_beta,
                             'a': alchemical_regions_pairs[0].softcore_a, 'b': alchemical_regions_pairs[0].softcore_b,
                             'c': alchemical_regions_pairs[0].softcore_c, 'd': alchemical_regions_pairs[0].softcore_d,
                             'e': alchemical_regions_pairs[0].softcore_e, 'f': alchemical_regions_pairs[0].softcore_f}

            # Get steric energy expressions.
            sterics_mixing_rules, exceptions_sterics_energy_expression = self._get_sterics_energy_expressions(lambda_var_suffixes, softcore_dict)

            # Define energy expression for sterics.
            sterics_energy_expression = exceptions_sterics_energy_expression + sterics_mixing_rules

            if not use_exact_pme_treatment:
                # There's no CustomNonbondedForce that models electrostatics if we use exact
                # PME treatment. Electrostatics is modeled through offset parameters.
                energy_expressions = self._get_electrostatics_energy_expressions(reference_force, lambda_var_suffixes, softcore_dict)
                (electrostatics_energy_expression,
                 exceptions_electrostatics_energy_expression) = energy_expressions  # Unpack tuple.

            # ------------------------------------------------------------
            # Create and configure all forces to add to alchemical system
            # ------------------------------------------------------------

            # Interactions and exceptions will be distributed according to the following table.

            # --------------------------------------------------------------------------------------------------
            # FORCE                                    | INTERACTION GROUP                                     |
            # --------------------------------------------------------------------------------------------------
            # nonbonded_force (unmodified)             | all interactions nonalchemical/nonalchemical          |
            #                                          | all exceptions nonalchemical/nonalchemical            |
            # --------------------------------------------------------------------------------------------------
            # aa_sterics_custom_nonbonded_force        | sterics interactions alchemical/alchemical            |
            # --------------------------------------------------------------------------------------------------
            # aa_electrostatics_custom_nonbonded_force | electrostatics interactions alchemical/alchemical     |
            #                                          | (only without exact PME treatment)                    |
            # --------------------------------------------------------------------------------------------------
            # na_sterics_custom_nonbonded_force        | sterics interactions non-alchemical/alchemical        |
            # --------------------------------------------------------------------------------------------------
            # na_electrostatics_custom_nonbonded_force | electrostatics interactions non-alchemical/alchemical |
            #                                          | (only without exact PME treatment)                    |
            # --------------------------------------------------------------------------------------------------
            # aa_sterics_custom_bond_force             | sterics exceptions alchemical/alchemical              |
            # --------------------------------------------------------------------------------------------------
            # aa_electrostatics_custom_bond_force      | electrostatics exceptions alchemical/alchemical       |
            #                                          | (only without exact PME treatment)                    |
            # --------------------------------------------------------------------------------------------------
            # na_sterics_custom_bond_force             | sterics exceptions non-alchemical/alchemical          |
            # --------------------------------------------------------------------------------------------------
            # na_electrostatics_custom_bond_force      | electrostatics exceptions non-alchemical/alchemical   |
            #                                          | (only without exact PME treatment)                    |
            # --------------------------------------------------------------------------------------------------

            def create_force(force_cls, energy_expression, lambda_variable_name, lambda_var_suffixes, is_lambda_controlled):
                """Shortcut to create a lambda-controlled custom forces."""
                if is_lambda_controlled:
                    force = force_cls(energy_expression)
                    for suffix in lambda_var_suffixes:
                        name = (lambda_variable_name + suffix)
                        force.addGlobalParameter(name, 1.0)
                else:  # fix lambda variable to 1.0
                    for suffix in lambda_var_suffixes:
                        name = (lambda_variable_name + suffix)
                        energy_expression = energy_expression + name + '=1.0;'
                    force = force_cls(energy_expression)
                return force

            # Create CustomNonbondedForces to handle sterics particle interactions between
            # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
            # to 1.0 for decoupled interactions in alchemical/alchemical force.
            if len(lambda_var_suffixes) > 1:
                aa_sterics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                                 'lambda_sterics', lambda_var_suffixes, is_lambda_controlled=True)
                all_sterics_custom_nonbonded_forces = [aa_sterics_custom_nonbonded_force]
            else:
                na_sterics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                                 'lambda_sterics', lambda_var_suffixes, is_lambda_controlled=True)
                aa_sterics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                                 'lambda_sterics', lambda_var_suffixes, alchemical_regions_pairs[0].annihilate_sterics)
                all_sterics_custom_nonbonded_forces = [na_sterics_custom_nonbonded_force, aa_sterics_custom_nonbonded_force]

            # Add parameters and configure CustomNonbondedForces to match reference force.
            for force in all_sterics_custom_nonbonded_forces:
                force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
                force.addPerParticleParameter("epsilon")  # Lennard-Jones epsilon
                force.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction())
                force.setCutoffDistance(nonbonded_force.getCutoffDistance())
                force.setSwitchingDistance(nonbonded_force.getSwitchingDistance())
                if self.disable_alchemical_dispersion_correction:
                    force.setUseLongRangeCorrection(False)
                else:
                    force.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())

                if is_periodic_method:
                    force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
                else:
                    force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

            if use_exact_pme_treatment:
                #electrostatics are handled by offset
                all_electrostatics_custom_nonbonded_forces = []
            else:
                # Create CustomNonbondedForces to handle electrostatics particle interactions between
                # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
                # to 1.0 for decoupled interactions in alchemical/alchemical force.
                if len(lambda_var_suffixes) > 1:
                    aa_electrostatics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                                            'lambda_electrostatics', lambda_var_suffixes, is_lambda_controlled=True)
                    all_electrostatics_custom_nonbonded_forces = [aa_electrostatics_custom_nonbonded_force]
                else:
                    na_electrostatics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                                            'lambda_electrostatics', lambda_var_suffixes, is_lambda_controlled=True)
                    aa_electrostatics_custom_nonbonded_force = create_force(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                                            'lambda_electrostatics', lambda_var_suffixes,  alchemical_regions_pairs[0].annihilate_electrostatics)
                    all_electrostatics_custom_nonbonded_forces = [na_electrostatics_custom_nonbonded_force,
                                                                  aa_electrostatics_custom_nonbonded_force]

            # Common parameters and configuration for electrostatics CustomNonbondedForces.
            for force in all_electrostatics_custom_nonbonded_forces:
                force.addPerParticleParameter("charge")  # partial charge
                force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
                if ((is_ewald_method and self.alchemical_pme_treatment == 'coulomb') or
                        (is_rf_method and self.alchemical_rf_treatment == 'switched')):
                    # Use switching function for alchemical electrostatics to ensure force continuity at cutoff.
                    force.setUseSwitchingFunction(True)
                else:
                    force.setUseSwitchingFunction(False)
                force.setSwitchingDistance(nonbonded_force.getCutoffDistance() - self.switch_width)
                force.setCutoffDistance(nonbonded_force.getCutoffDistance())
                force.setUseLongRangeCorrection(False)  # long-range dispersion correction is meaningless for electrostatics

                if is_periodic_method:
                    force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
                else:
                    force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

            # Create CustomBondForces to handle sterics 1,4 exceptions interactions between
            # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
            # to 1.0 for decoupled interactions in alchemical/alchemical force.
            if len(lambda_var_suffixes) > 1:
                aa_sterics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                            'lambda_sterics', lambda_var_suffixes, is_lambda_controlled=True)
                all_sterics_custom_bond_forces = [aa_sterics_custom_bond_force]
            else:
                na_sterics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                            'lambda_sterics', lambda_var_suffixes, is_lambda_controlled=True)
                aa_sterics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                            'lambda_sterics', lambda_var_suffixes, alchemical_regions_pairs[0].annihilate_sterics)
                all_sterics_custom_bond_forces = [na_sterics_custom_bond_force, aa_sterics_custom_bond_force]

            for force in all_sterics_custom_bond_forces:
                force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma
                force.addPerBondParameter("epsilon")  # Lennard-Jones effective epsilon

            # With exact PME treatment, exception electrostatics is handled through offset parameters.
            if use_exact_pme_treatment:
                all_electrostatics_custom_bond_forces = []
            else:
                # Create CustomBondForces to handle electrostatics 1,4 exceptions interactions between
                # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
                # to 1.0 for decoupled interactions in alchemical/alchemical force.
                if len(lambda_var_suffixes) > 1:
                    aa_electrostatics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                                       'lambda_electrostatics', lambda_var_suffixes, is_lambda_controlled=True)
                    all_electrostatics_custom_bond_forces = [aa_electrostatics_custom_bond_force]
                else:
                    na_electrostatics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                                       'lambda_electrostatics', lambda_var_suffixes, is_lambda_controlled=True)
                    aa_electrostatics_custom_bond_force = create_force(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                                       'lambda_electrostatics', lambda_var_suffixes, alchemical_regions_pairs[0].annihilate_electrostatics)
                    all_electrostatics_custom_bond_forces = [na_electrostatics_custom_bond_force, aa_electrostatics_custom_bond_force]

            # Create CustomBondForce to handle exceptions for electrostatics
            for force in all_electrostatics_custom_bond_forces:
                force.addPerBondParameter("chargeprod")  # charge product
                force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma

            # -------------------------------------------------------------------------------
            # Distribute particle interactions contributions in appropriate nonbonded forces
            # -------------------------------------------------------------------------------

            # Create atom groups.
            alchemical_atomsets = [region.alchemical_atoms for region in alchemical_regions_pairs]

            # Copy NonbondedForce particle terms for alchemically-modified particles
            # to CustomNonbondedForces, and/or add the charge offsets for exact PME.
            # On CUDA, for efficiency reasons, all nonbonded forces (custom and not)
            # must have the same particles.
            for particle_index in range(nonbonded_force.getNumParticles()):
                # Retrieve nonbonded parameters.
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
                # Set sterics parameters in CustomNonbondedForces.
                for force in all_sterics_custom_nonbonded_forces:
                    force.addParticle([sigma, epsilon])
                # Set electrostatics parameters in CustomNonbondedForces.
                for force in all_electrostatics_custom_nonbonded_forces:
                    force.addParticle([charge, sigma])
                # Set offset parameters in NonbondedForce.
                if use_exact_pme_treatment and particle_index in alchemical_atomsets[0] and len(lambda_var_suffixes) == 1:
                    nonbonded_force.addParticleParameterOffset('lambda_electrostatics{}'.format(lambda_var_suffixes[0]),
                                                               particle_index, charge, 0.0, 0.0)

            # Turn off interactions contribution from alchemically-modified particles in unmodified
            # NonbondedForce that will be handled by all other forces
            for particle_index in range(nonbonded_force.getNumParticles()):
                # Retrieve parameters.
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
                # Even with exact treatment of the PME electrostatics, we turn off
                # the NonbondedForce charge which is modeled by the offset parameter.
                if particle_index in alchemical_atomsets[0]:
                    nonbonded_force.setParticleParameters(particle_index, abs(0.0*charge), sigma, abs(0*epsilon))

            # Restrict interaction evaluation of CustomNonbondedForces to their respective atom groups.
            # Sterics
            if len(lambda_var_suffixes) == 1:
                logger.debug('Adding steric interaction groups between {} and the environment.'.format(lambda_var_suffixes[0]))
                na_sterics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomsets[0])

            logger.debug('Adding a steric interaction group between group {0} and {1}.'.format(lambda_var_suffixes[0],
                                                                                               lambda_var_suffixes[-1]))
            aa_sterics_custom_nonbonded_force.addInteractionGroup(alchemical_atomsets[0], alchemical_atomsets[-1])

            # Electrostatics
            if not use_exact_pme_treatment:
                if len(lambda_var_suffixes) == 1:
                    logger.debug('Adding electrostatic interaction groups between {} and the environment.'.format(lambda_var_suffixes[0]))
                    na_electrostatics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomsets[0])

                logger.debug('Adding a electrostatic interaction group between group {0} and {1}.'.format(lambda_var_suffixes[0],
                                                                                                          lambda_var_suffixes[-1]))
                aa_electrostatics_custom_nonbonded_force.addInteractionGroup(alchemical_atomsets[0], alchemical_atomsets[-1])

            else:
                # Using the nonbonded force to handle electrostatics
                # and the "interaction groups" in the nonbonded have already been handled by exclusions.
                pass

            # ---------------------------------------------------------------
            # Distribute exceptions contributions in appropriate bond forces
            # ---------------------------------------------------------------

            all_custom_nonbonded_forces = all_sterics_custom_nonbonded_forces + all_electrostatics_custom_nonbonded_forces

            # Move all NonbondedForce exception terms for alchemically-modified particles to CustomBondForces.
            for exception_index in range(nonbonded_force.getNumExceptions()):
                # Retrieve parameters.
                iatom, jatom, chargeprod, sigma, epsilon = nonbonded_force.getExceptionParameters(exception_index)

                # Exclude this atom pair in CustomNonbondedForces. All nonbonded forces
                # must have the same number of exceptions/exclusions on CUDA platform.
                for force in all_custom_nonbonded_forces:
                    force.addExclusion(iatom, jatom)

                # Check if this is an exception or an exclusion
                is_exception_epsilon = abs(epsilon.value_in_unit_system(unit.md_unit_system)) > 0.0
                is_exception_chargeprod = abs(chargeprod.value_in_unit_system(unit.md_unit_system)) > 0.0

                # Check how many alchemical atoms we have in the exception.
                if len(lambda_var_suffixes) > 1:
                    # Pair of interacting alchemical regions, therefore they are both alchemical or neither alchemical.
                    both_alchemical = ((iatom in alchemical_atomsets[0] and jatom in alchemical_atomsets[1]) or
                                       (jatom in alchemical_atomsets[0] and iatom in alchemical_atomsets[1]))
                    only_one_alchemical = False
                    #The condition of at_least_one_alchemical is treated only once per single
                    # region so we don't repeat it when dealing with pairs of interacting regions.
                    at_least_one_alchemical = False


                    if use_exact_pme_treatment and both_alchemical and is_exception_chargeprod:
                        # Exceptions here should be scaled by lam0*lam1.
                        # This can be implemented in the future using a CustomBondForce.
                        raise ValueError('Cannot have exception that straddles two alchemical regions')

                else:
                    # Single alchemical region.
                    both_alchemical = iatom in alchemical_atomsets[0] and jatom in alchemical_atomsets[0]
                    at_least_one_alchemical = iatom in alchemical_atomsets[0] or jatom in alchemical_atomsets[0]
                    only_one_alchemical = at_least_one_alchemical and not both_alchemical

                    # If this is an electrostatic exception and we're using exact PME,
                    # we just have to add the exception offset to the NonbondedForce.
                    if use_exact_pme_treatment and at_least_one_alchemical and is_exception_chargeprod:
                        nonbonded_force.addExceptionParameterOffset('lambda_electrostatics{}'.format(lambda_var_suffixes[0]),
                                                                    exception_index, chargeprod, 0.0, 0.0)

                # If exception (and not exclusion), add special CustomBondForce terms to
                # handle alchemically-modified Lennard-Jones and electrostatics exceptions
                if both_alchemical:
                    if is_exception_epsilon:
                        aa_sterics_custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])
                    if is_exception_chargeprod and not use_exact_pme_treatment:
                        aa_electrostatics_custom_bond_force.addBond(iatom, jatom, [chargeprod, sigma])

                # When this is a single region we model the exception between alchemical
                # and non-alchemical particles using a single custom bond.
                elif only_one_alchemical:
                    if is_exception_epsilon:
                        na_sterics_custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])
                    if is_exception_chargeprod and not use_exact_pme_treatment:
                        na_electrostatics_custom_bond_force.addBond(iatom, jatom, [chargeprod, sigma])
                # else: both particles are non-alchemical, leave them in the unmodified NonbondedForce

                # Turn off all exception contributions from alchemical atoms in the NonbondedForce
                # modelling non-alchemical atoms only. We need to do it only once per single
                # region so we don't repeat it when dealing with pairs of interacting regions.
                if at_least_one_alchemical:
                    nonbonded_force.setExceptionParameters(exception_index, iatom, jatom,
                                                           abs(0.0*chargeprod), sigma, abs(0.0*epsilon))

            # With exact treatment of PME electrostatics, the NonbondedForce
            # is affected by lambda electrostatics as well.
            if 'lambda_sterics{}'.format(lambda_var_suffixes[0]) in forces_by_lambda:
                forces_by_lambda['lambda_electrostatics{}'.format(lambda_var_suffixes[0])].extend(all_electrostatics_custom_nonbonded_forces + all_electrostatics_custom_bond_forces)
                forces_by_lambda['lambda_sterics{}'.format(lambda_var_suffixes[0])].extend(all_sterics_custom_nonbonded_forces + all_sterics_custom_bond_forces)
            else:
                forces_by_lambda['lambda_electrostatics{}'.format(lambda_var_suffixes[0])] = all_electrostatics_custom_nonbonded_forces + all_electrostatics_custom_bond_forces
                forces_by_lambda['lambda_sterics{}'.format(lambda_var_suffixes[0])] = all_sterics_custom_nonbonded_forces + all_sterics_custom_bond_forces

        if use_exact_pme_treatment:
            forces_by_lambda['lambda_electrostatics{}'.format(lambda_var_suffixes[0])].append(nonbonded_force)
        else:
            forces_by_lambda[''] = [nonbonded_force]
        return forces_by_lambda
