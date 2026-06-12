import cobra.io as io
from optlang import Constraint, Model, Objective, Variable


class Simulator:
    """
    Constraint-based metabolic model simulator implementing FBA, Linear MOMA, and ROOM methods.

    This class provides methods for:
    - Flux Balance Analysis (FBA)
    - Linear Minimization of Metabolic Adjustment (Linear MOMA)
    - Regulatory On/Off Minimization (ROOM)
    """

    def __init__(self):
        """
        Constructor for Simulator

        Initializes all model components to None. These will be populated
        when a model is loaded using read_model() or load_cobra_model().
        """
        self.cobra_model = None
        self.model_metabolites = None
        self.model_reactions = None
        self.model_genes = None  # Added to avoid AttributeError
        self.Smatrix = None
        self.lower_boundary_constraints = None
        self.upper_boundary_constraints = None
        self.objective = None

    def run_MOMA(self, wild_flux={}, flux_constraints={}, inf_flag=False):
        """
        Linear MOMA (Minimization of Metabolic Adjustment) analysis.

        Linear MOMA finds a flux distribution that minimizes the Manhattan distance (L1 norm)
        to a reference (wild-type) flux distribution while satisfying stoichiometric
        and flux bound constraints. This is computationally more efficient than
        quadratic MOMA while providing similar biological insights.

        Parameters:
        -----------
        wild_flux : dict
            Reference flux distribution (typically wild-type). Keys are reaction IDs,
            values are flux values.
        flux_constraints : dict
            Additional flux constraints. Keys are reaction IDs, values are tuples
            of (lower_bound, upper_bound).
        inf_flag : bool, optional (default=False)
            If False, replaces infinite bounds with ±1000.

        Returns:
        --------
        tuple: (status, objective_value, flux_distribution)
            - status: Optimization status ('optimal' if successful)
            - objective_value: Minimized Manhattan distance (sum of absolute deviations)
            - flux_distribution: Dict of reaction IDs to flux values

        Notes:
        ------
        The objective minimizes: sum(|v - w|) where v is mutant flux and w is wild-type flux.
        This is implemented using auxiliary variables: |v - w| = fplus + fminus
        where v - w = fplus - fminus, fplus >= 0, fminus >= 0.
        """
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        Smatrix = self.Smatrix
        lower_boundary_constraints = self.lower_boundary_constraints.copy()
        upper_boundary_constraints = self.upper_boundary_constraints.copy()

        # Replace infinite bounds with finite values if inf_flag is False
        if not inf_flag:
            for key in lower_boundary_constraints:
                if lower_boundary_constraints[key] == float("-inf"):
                    lower_boundary_constraints[key] = -1000.0

            for key in upper_boundary_constraints:
                if upper_boundary_constraints[key] == float("inf"):
                    upper_boundary_constraints[key] = 1000.0

        # Create optlang model
        m = Model(name="Linear_MOMA")

        # Create variables
        v = {}  # Flux variables
        fplus = {}  # Positive deviation from wild-type
        fminus = {}  # Negative deviation from wild-type

        variables_to_add = []

        for each_reaction in model_reactions:
            # Set flux bounds based on constraints or model defaults
            if each_reaction in flux_constraints:
                v[each_reaction] = Variable(
                    each_reaction, lb=flux_constraints[each_reaction][0], ub=flux_constraints[each_reaction][1]
                )
            else:
                v[each_reaction] = Variable(
                    each_reaction,
                    lb=lower_boundary_constraints[each_reaction],
                    ub=upper_boundary_constraints[each_reaction],
                )
            # Deviation variables (non-negative)
            fplus[each_reaction] = Variable(f"fplus_{each_reaction}", lb=0.0, ub=1000.0)
            fminus[each_reaction] = Variable(f"fminus_{each_reaction}", lb=0.0, ub=1000.0)

            variables_to_add.extend([v[each_reaction], fplus[each_reaction], fminus[each_reaction]])

        m.add(variables_to_add)

        # Add constraints for absolute value decomposition
        # For each reaction: v - w = fplus - fminus
        # This means: |v - w| = fplus + fminus (when both are minimized)
        constraints = []

        for each_reaction in model_reactions:
            if each_reaction in wild_flux:
                w = wild_flux[each_reaction]  # Wild-type flux

                # Constraint: v - w = fplus - fminus
                # Rearranged: v - fplus + fminus = w
                constraints.append(
                    Constraint(
                        v[each_reaction] - fplus[each_reaction] + fminus[each_reaction],
                        lb=w,
                        ub=w,
                        name=f"abs_decomp_{each_reaction}",
                    )
                )

        m.add(constraints)

        # Add steady-state mass balance constraints (Sv = 0)
        mass_balance_constraints = []

        for each_metabolite in model_metabolites:
            # Find all reactions involving this metabolite
            metabolite_reactions = [(met, rxn) for (met, rxn) in Smatrix.keys() if met == each_metabolite]

            if len(metabolite_reactions) == 0:
                continue

            # Create mass balance expression: sum(S_ij * v_j) = 0
            expr = sum(v[reaction] * Smatrix[metabolite, reaction] for metabolite, reaction in metabolite_reactions)

            mass_balance_constraints.append(Constraint(expr, lb=0, ub=0, name=f"mass_balance_{each_metabolite}"))

        m.add(mass_balance_constraints)

        # Set objective: minimize sum of absolute deviations (Manhattan distance)
        target_reactions = wild_flux.keys()
        objective_expr = sum((fplus[each_reaction] + fminus[each_reaction]) for each_reaction in target_reactions)
        m.objective = Objective(objective_expr, direction="min")

        # Solve the optimization problem
        m.optimize()

        # Extract results if optimization was successful
        if m.status == "optimal":
            flux_distribution = {}
            for reaction in model_reactions:
                flux_distribution[reaction] = float(v[reaction].primal)
                # Round near-zero values to exactly zero
                if abs(float(v[reaction].primal)) <= 1e-6:
                    flux_distribution[reaction] = 0.0

            return m.status, m.objective.value, flux_distribution
        else:
            return m.status, False, False

    def run_ROOM(self, wild_flux={}, flux_constraints={}, delta=0.03, epsilon=0.001, inf_flag=False):
        """
        Regulatory On/Off Minimization (ROOM) analysis.

        ROOM finds a flux distribution that minimizes the number of significant flux changes
        compared to a reference (wild-type) flux distribution. It uses binary variables to
        count reactions that deviate significantly from the reference state.

        Parameters:
        -----------
        wild_flux : dict
            Reference flux distribution (typically wild-type). Keys are reaction IDs,
            values are flux values.
        flux_constraints : dict
            Additional flux constraints. Keys are reaction IDs, values are tuples
            of (lower_bound, upper_bound).
        delta : float, optional (default=0.03)
            Relative tolerance for flux changes. A reaction is considered "changed"
            if |v - w| > delta * |w| where v is mutant flux and w is wild-type flux.
        epsilon : float, optional (default=0.001)
            Absolute tolerance for near-zero fluxes. Used to avoid division by zero
            and to define inactive reactions.
        inf_flag : bool, optional (default=False)
            If False, replaces infinite bounds with ±1000.

        Returns:
        --------
        tuple: (status, objective_value, flux_distribution)
            - status: Gurobi optimization status (2 = optimal)
            - objective_value: Number of significantly changed reactions
            - flux_distribution: Dict of reaction IDs to flux values

        Notes:
        ------
        ROOM is particularly useful for predicting gene knockout phenotypes, as it
        assumes the cell minimizes regulatory changes rather than metabolic adjustment.
        """
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        Smatrix = self.Smatrix
        lower_boundary_constraints = self.lower_boundary_constraints.copy()
        upper_boundary_constraints = self.upper_boundary_constraints.copy()

        # Replace infinite bounds with finite values if inf_flag is False
        if not inf_flag:
            for key in lower_boundary_constraints:
                if lower_boundary_constraints[key] == float("-inf"):
                    lower_boundary_constraints[key] = -1000.0

            for key in upper_boundary_constraints:
                if upper_boundary_constraints[key] == float("inf"):
                    upper_boundary_constraints[key] = 1000.0

        # Create optlang model
        m = Model(name="ROOM")

        # Create variables
        v = {}  # Flux variables
        y = {}  # Binary variables: 1 if reaction flux significantly changed, 0 otherwise

        variables_to_add = []

        for each_reaction in model_reactions:
            # Set flux bounds based on constraints or model defaults
            if each_reaction in flux_constraints:
                v[each_reaction] = Variable(
                    each_reaction, lb=flux_constraints[each_reaction][0], ub=flux_constraints[each_reaction][1]
                )
            else:
                v[each_reaction] = Variable(
                    each_reaction,
                    lb=lower_boundary_constraints[each_reaction],
                    ub=upper_boundary_constraints[each_reaction],
                )

            # Binary indicator variable for significant flux change
            y[each_reaction] = Variable(f"y_{each_reaction}", type="binary")

            variables_to_add.extend([v[each_reaction], y[each_reaction]])

        m.add(variables_to_add)

        # Add steady-state mass balance constraints (Sv = 0)
        mass_balance_constraints = []

        for each_metabolite in model_metabolites:
            # Find all reactions involving this metabolite
            metabolite_reactions = [(met, rxn) for (met, rxn) in Smatrix.keys() if met == each_metabolite]

            if len(metabolite_reactions) == 0:
                continue

            # Create mass balance expression: sum(S_ij * v_j) = 0
            expr = sum(v[reaction] * Smatrix[metabolite, reaction] for metabolite, reaction in metabolite_reactions)

            mass_balance_constraints.append(Constraint(expr, lb=0, ub=0, name=f"mass_balance_{each_metabolite}"))

        m.add(mass_balance_constraints)

        # Add ROOM-specific constraints
        # For each reaction with reference flux, constrain deviation based on binary variable
        for each_reaction in model_reactions:
            if each_reaction not in wild_flux:
                continue

            w = wild_flux[each_reaction]  # Reference (wild-type) flux

            # Calculate upper and lower bounds for allowable flux range
            if abs(w) < epsilon:
                # If wild-type flux is near zero, use absolute tolerance
                w_upper = epsilon
                w_lower = -epsilon
            else:
                # Use relative tolerance based on wild-type flux magnitude
                w_upper = w + delta * abs(w)
                w_lower = w - delta * abs(w)

            # Big-M formulation to link binary variable y to flux deviation
            # If y = 0 (no significant change), then w_lower <= v <= w_upper
            # If y = 1 (significant change allowed), then lb <= v <= ub (no additional constraint)

            # Get actual bounds for the reaction
            if each_reaction in flux_constraints:
                lb = flux_constraints[each_reaction][0]
                ub = flux_constraints[each_reaction][1]
            else:
                lb = lower_boundary_constraints[each_reaction]
                ub = upper_boundary_constraints[each_reaction]

            # Big-M values (should be larger than any feasible flux)
            M_upper = ub - w_upper if ub < 1000 else 2000
            M_lower = w_lower - lb if lb > -1000 else 2000

            # v <= w_upper + M_upper * y
            # When y=0: v <= w_upper (enforces upper bound of allowable range)
            # When y=1: v <= w_upper + M_upper (effectively no constraint)
            m.add(
                Constraint(
                    v[each_reaction] - M_upper * y[each_reaction], ub=w_upper, name=f"room_upper_{each_reaction}"
                )
            )

            # v >= w_lower - M_lower * y
            # When y=0: v >= w_lower (enforces lower bound of allowable range)
            # When y=1: v >= w_lower - M_lower (effectively no constraint)
            m.add(
                Constraint(
                    v[each_reaction] + M_lower * y[each_reaction], lb=w_lower, name=f"room_lower_{each_reaction}"
                )
            )

        # Set objective: minimize number of significantly changed reactions
        objective_expr = sum(y[each_reaction] for each_reaction in wild_flux.keys())
        m.objective = Objective(objective_expr, direction="min")

        # Solve the optimization problem
        m.optimize()

        # Extract results if optimization was successful
        if m.status == "optimal":
            objective_value = m.objective.value  # Number of changed reactions
            flux_distribution = {}
            for reaction in model_reactions:
                flux_distribution[reaction] = float(v[reaction].primal)
                # Round near-zero values to exactly zero
                if abs(float(v[reaction].primal)) <= 1e-6:
                    flux_distribution[reaction] = 0.0

            return m.status, objective_value, flux_distribution
        else:
            return m.status, False, False

    def run_FSEOF(self, target_reaction, objective_fraction=1.0, flux_constraints={}, n_points=10, inf_flag=False):
        """
        Flux Scanning based on Enforced Objective Flux (FSEOF).

        Parameters:
        -----------
        target_reaction : str
            Reaction ID to scan
        objective_fraction : float
            Fraction of maximum objective flux to maintain
        n_points : int
            Number of scanning points between min and max flux
        flux_constraints : dict
            Additional flux constraints
        inf_flag : bool
            If False, replaces infinite bounds with ±1000

        Returns:
        --------
        tuple: (status, results_df, correlations_df)
            - status: 'optimal' or 'failed'
            - results_df: DataFrame with columns [target_flux, objective_flux, reaction1, reaction2, ...]
            - correlations_df: DataFrame with correlation coefficients
        """
        import numpy as np
        import pandas as pd

        # Step 1: Find maximum objective value
        fba_status, max_obj_value, _ = self.run_FBA(flux_constraints=flux_constraints, inf_flag=inf_flag)

        if fba_status != "optimal":
            return fba_status, None, None

        # Step 2: Set minimum objective constraint
        min_obj_flux = objective_fraction * max_obj_value
        objective = self.objective
        fseof_constraints = flux_constraints.copy()
        fseof_constraints[objective] = (min_obj_flux, float("inf"))

        # Step 3: Find target reaction's feasible range
        min_status, min_flux_val, _ = self.run_FBA(
            new_objective=target_reaction, flux_constraints=fseof_constraints, inf_flag=inf_flag, mode="min"
        )

        max_status, max_flux_val, _ = self.run_FBA(
            new_objective=target_reaction, flux_constraints=fseof_constraints, inf_flag=inf_flag, mode="max"
        )

        if min_status != "optimal" or max_status != "optimal":
            return "failed", None, None

        # Step 4: Generate scan points
        scan_points = np.linspace(min_flux_val, max_flux_val, n_points)

        # Step 5: Scan through each point
        flux_data = []

        for target_flux in scan_points:
            # Fix target reaction at current scan point
            scan_constraints = fseof_constraints.copy()
            scan_constraints[target_reaction] = (target_flux, target_flux)

            # Run FBA
            status, obj_val, flux_dist = self.run_FBA(flux_constraints=scan_constraints, inf_flag=inf_flag)

            if status == "optimal":
                # Store flux distribution with target and objective values
                row_data = flux_dist.copy()
                row_data["target_flux"] = target_flux
                row_data["objective_flux"] = obj_val
                flux_data.append(row_data)

        if not flux_data:
            return "failed", None, None

        # Step 6: Create results DataFrame
        results_df = pd.DataFrame(flux_data)

        # Reorder columns: target_flux, objective_flux, then all reactions
        cols = ["target_flux", "objective_flux"] + [
            col for col in results_df.columns if col not in ["target_flux", "objective_flux"]
        ]
        results_df = results_df[cols]

        # Step 7: Calculate correlations with target reaction
        reaction_columns = [col for col in results_df.columns if col not in ["target_flux", "objective_flux"]]

        correlations = {}
        for rxn in reaction_columns:
            corr = results_df["target_flux"].corr(results_df[rxn])
            correlations[rxn] = corr

        # Create correlations DataFrame
        correlations_df = (
            pd.DataFrame(
                [
                    {"reaction": rxn, "correlation": corr, "abs_correlation": abs(corr)}
                    for rxn, corr in correlations.items()
                ]
            )
            .sort_values("abs_correlation", ascending=False)
            .reset_index(drop=True)
        )

        return "optimal", results_df, correlations_df

    def run_FVA(self, fraction_of_optimum=1.0, flux_constraints={}, inf_flag=False, reactions_to_analyze=None):
        """
        Flux Variability Analysis (FVA).

        FVA determines the minimum and maximum possible flux through each reaction
        while maintaining a minimum objective flux. This reveals the flexibility
        of the metabolic network.

        Parameters:
        -----------
        fraction_of_optimum : float, optional (default=1.0)
            Fraction of maximum objective flux to maintain (0 to 1).
        flux_constraints : dict
            Additional flux constraints. Keys are reaction IDs, values are tuples
            of (lower_bound, upper_bound).
        inf_flag : bool, optional (default=False)
            If False, replaces infinite boundobjective_fractions with ±1000.
        reactions_to_analyze : list, optional
            List of reaction IDs to analyze. If None, analyzes all reactions.

        Returns:
        --------
        dict: FVA results containing:
            - 'status': 'optimal' if successful
            - 'fva_data': Dict mapping reaction IDs to {'minimum': min_flux, 'maximum': max_flux}
            - 'objective_value': The enforced objective value
        """
        # First, find maximum objective value
        fba_status, max_obj_value, _ = self.run_FBA(flux_constraints=flux_constraints, inf_flag=inf_flag)

        if fba_status != "optimal":
            return {"status": fba_status, "fva_data": {}, "objective_value": 0}

        # Calculate minimum required objective flux
        min_obj_flux = fraction_of_optimum * max_obj_value

        # Determine which reactions to analyze
        if reactions_to_analyze is None:
            reactions_to_analyze = self.model_reactions

        # Add constraint to maintain minimum objective
        objective = self.objective
        fva_constraints = flux_constraints.copy()
        fva_constraints[objective] = (min_obj_flux, float("inf"))

        # Perform FVA for each reaction
        fva_data = {}

        for reaction in reactions_to_analyze:
            # Minimize reaction flux
            min_status, min_flux, _ = self.run_FBA(
                new_objective=reaction, flux_constraints=fva_constraints, inf_flag=inf_flag, mode="min"
            )

            # Maximize reaction flux
            max_status, max_flux, _ = self.run_FBA(
                new_objective=reaction, flux_constraints=fva_constraints, inf_flag=inf_flag, mode="max"
            )

            if min_status == "optimal" and max_status == "optimal":
                fva_data[reaction] = {"minimum": min_flux, "maximum": max_flux, "range": max_flux - min_flux}

        return {"status": "optimal", "fva_data": fva_data, "objective_value": min_obj_flux}

    def run_FBA(
        self, new_objective="", flux_constraints={}, inf_flag=False, internal_flux_minimization=False, mode="max"
    ):
        """
        Flux Balance Analysis (FBA).

        FBA maximizes or minimizes an objective function (typically biomass or a target
        metabolite production) subject to stoichiometric constraints and flux bounds.

        Parameters:
        -----------
        new_objective : str, optional
            Reaction ID to use as objective function. If empty, uses model's default objective.
        flux_constraints : dict
            Additional flux constraints. Keys are reaction IDs, values are tuples
            of (lower_bound, upper_bound).
        inf_flag : bool, optional (default=False)
            If False, replaces infinite bounds with ±1000.
        internal_flux_minimization : bool, optional (default=False)
            If True, performs parsimonious FBA (pFBA) to minimize total flux while
            maintaining optimal objective value.
        mode : str, optional (default='max')
            Optimization mode: 'max' for maximization, 'min' for minimization.

        Returns:
        --------
        tuple: (status, objective_value, flux_distribution)
            - status: Gurobi optimization status (2 = optimal)
            - objective_value: Optimized objective function value (or total flux for pFBA)
            - flux_distribution: Dict of reaction IDs to flux values

        Notes:
        ------
        Parsimonious FBA (pFBA) is useful for obtaining more biologically realistic
        flux distributions by minimizing the total sum of fluxes while maintaining
        optimal growth or production.
        """
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        # Determine objective reaction
        if new_objective == "":
            objective = self.objective
        else:
            objective = new_objective

        Smatrix = self.Smatrix
        lower_boundary_constraints = self.lower_boundary_constraints.copy()
        upper_boundary_constraints = self.upper_boundary_constraints.copy()

        # Replace infinite bounds with finite values if inf_flag is False
        if not inf_flag:
            for key in lower_boundary_constraints:
                if lower_boundary_constraints[key] == float("-inf"):
                    lower_boundary_constraints[key] = -1000.0

            for key in upper_boundary_constraints:
                if upper_boundary_constraints[key] == float("inf"):
                    upper_boundary_constraints[key] = 1000.0

        # Create optlang model
        m = Model(name="FBA")

        # Create variables
        v = {}  # Flux variables
        fplus = {}  # Positive flux (for pFBA)
        fminus = {}  # Negative flux (for pFBA)

        variables_to_add = []

        for each_reaction in model_reactions:
            # Set flux bounds based on constraints or model defaults
            if each_reaction in flux_constraints:
                v[each_reaction] = Variable(
                    each_reaction, lb=flux_constraints[each_reaction][0], ub=flux_constraints[each_reaction][1]
                )
            else:
                v[each_reaction] = Variable(
                    each_reaction,
                    lb=lower_boundary_constraints[each_reaction],
                    ub=upper_boundary_constraints[each_reaction],
                )
            # Auxiliary variables for pFBA
            fplus[each_reaction] = Variable(f"fplus_{each_reaction}", lb=0.0, ub=1000.0)
            fminus[each_reaction] = Variable(f"fminus_{each_reaction}", lb=0.0, ub=1000.0)

            variables_to_add.extend([v[each_reaction], fplus[each_reaction], fminus[each_reaction]])

        m.add(variables_to_add)

        # Add steady-state mass balance constraints (Sv = 0)
        mass_balance_constraints = []

        for each_metabolite in model_metabolites:
            # Find all reactions involving this metabolite
            metabolite_reactions = [(met, rxn) for (met, rxn) in Smatrix.keys() if met == each_metabolite]

            if len(metabolite_reactions) == 0:
                continue

            # Create mass balance expression: sum(S_ij * v_j) = 0
            expr = sum(v[reaction] * Smatrix[metabolite, reaction] for metabolite, reaction in metabolite_reactions)

            mass_balance_constraints.append(Constraint(expr, lb=0, ub=0, name=f"mass_balance_{each_metabolite}"))

        m.add(mass_balance_constraints)

        # Set primary objective (maximize or minimize target reaction flux)
        if mode == "max":
            m.objective = Objective(v[objective], direction="max")
        elif mode == "min":
            m.objective = Objective(v[objective], direction="min")

        # Solve primary optimization
        m.optimize()

        if m.status == "optimal":
            objective_value = m.objective.value

            # Parsimonious FBA: minimize total flux while maintaining optimal objective
            if internal_flux_minimization:
                # Fix objective flux to optimal value
                m.add(Constraint(v[objective], lb=objective_value, ub=objective_value, name="fix_objective"))

                # Add flux decomposition constraints for all reactions
                pfba_constraints = []
                for each_reaction in model_reactions:
                    pfba_constraints.append(
                        Constraint(
                            fplus[each_reaction] - fminus[each_reaction] - v[each_reaction],
                            lb=0,
                            ub=0,
                            name=f"flux_decomp_{each_reaction}",
                        )
                    )

                m.add(pfba_constraints)

                # New objective: minimize sum of absolute fluxes
                objective_expr = sum(
                    (fplus[each_reaction] + fminus[each_reaction]) for each_reaction in model_reactions
                )
                m.objective = Objective(objective_expr, direction="min")

                # Solve secondary optimization
                m.optimize()

                if m.status == "optimal":
                    objective_value = m.objective.value  # Total flux sum
                    flux_distribution = {}
                    for reaction in model_reactions:
                        flux_distribution[reaction] = float(v[reaction].primal)
                        # Round near-zero values to exactly zero
                        if abs(float(v[reaction].primal)) <= 1e-6:
                            flux_distribution[reaction] = 0.0
                    return m.status, objective_value, flux_distribution
                else:
                    return m.status, False, False
            else:
                # Standard FBA: return optimal flux distribution
                flux_distribution = {}
                for reaction in model_reactions:
                    flux_distribution[reaction] = float(v[reaction].primal)
                    # Round near-zero values to exactly zero
                    if abs(float(v[reaction].primal)) <= 1e-6:
                        flux_distribution[reaction] = 0.0
                return m.status, objective_value, flux_distribution

        return m.status, False, False

    def read_model(self, filename):
        """
        Load a metabolic model from an SBML file.

        Parameters:
        -----------
        filename : str
            Path to the SBML model file.

        Returns:
        --------
        tuple: (metabolites, reactions, Smatrix, lb, ub, objective)
            Components of the loaded model.
        """
        model = io.read_sbml_model(filename)
        return self.load_cobra_model(model)

    def load_cobra_model(self, cobra_model):
        """
        Load a COBRApy model into the simulator.

        Extracts all necessary components from a COBRApy model including metabolites,
        reactions, stoichiometric matrix, flux bounds, and objective function.

        Parameters:
        -----------
        cobra_model : cobra.Model
            A COBRApy model object.

        Returns:
        --------
        tuple: (metabolites, reactions, Smatrix, lb, ub, objective)
            - metabolites: List of metabolite IDs
            - reactions: List of reaction IDs
            - Smatrix: Dict mapping (metabolite, reaction) to stoichiometric coefficient
            - lb: Dict of lower bounds for each reaction
            - ub: Dict of upper bounds for each reaction
            - objective: Reaction ID of the objective function
        """
        self.cobra_model = cobra_model
        model = cobra_model
        model_metabolites = []
        model_reactions = []
        model_genes = []
        lower_boundary_constraints = {}
        upper_boundary_constraints = {}
        objective_reaction = ""

        # Extract metabolites
        for each_metabolite in model.metabolites:
            model_metabolites.append(each_metabolite.id)

        # Extract genes
        model_genes = [each_gene.id for each_gene in model.genes]

        # Build stoichiometric matrix
        Smatrix = {}

        for each_reaction in model.reactions:
            # Identify objective reaction
            if each_reaction.objective_coefficient == 1.0:
                objective_reaction = each_reaction.id

            # Extract stoichiometric coefficients for reactants
            reactant_list = each_reaction.reactants
            reactant_coff_list = each_reaction.get_coefficients(reactant_list)

            # Extract stoichiometric coefficients for products
            product_list = each_reaction.products
            product_coff_list = each_reaction.get_coefficients(product_list)

            reactant_coff_list = list(reactant_coff_list)
            product_coff_list = list(product_coff_list)

            # Add reactant coefficients to Smatrix (negative by convention)
            for i in range(len(reactant_list)):
                Smatrix[(reactant_list[i].id, each_reaction.id)] = reactant_coff_list[i]

            # Add product coefficients to Smatrix (positive by convention)
            for i in range(len(product_list)):
                Smatrix[(product_list[i].id, each_reaction.id)] = product_coff_list[i]

            # Store reaction ID
            model_reactions.append(each_reaction.id)

            # Extract flux bounds
            lb = each_reaction.lower_bound
            ub = each_reaction.upper_bound

            # Convert very large bounds to infinity for cleaner representation
            if lb < -1000.0:
                lb = float("-inf")
            if ub > 1000.0:
                ub = float("inf")

            lower_boundary_constraints[each_reaction.id] = lb
            upper_boundary_constraints[each_reaction.id] = ub

        # Store all model components as instance variables
        self.model_metabolites = model_metabolites
        self.model_reactions = model_reactions
        self.model_genes = model_genes
        self.Smatrix = Smatrix
        self.lower_boundary_constraints = lower_boundary_constraints
        self.upper_boundary_constraints = upper_boundary_constraints
        self.objective = objective_reaction

        return (
            model_metabolites,
            model_reactions,
            Smatrix,
            lower_boundary_constraints,
            upper_boundary_constraints,
            objective_reaction,
        )


def test_simulator_comparison(model_name="iJO1366"):
    """
    Test function to compare Simulator (optlang) results with COBRApy results
    using E. coli model for FBA, Linear MOMA, and ROOM analyses.

    Parameters:
    -----------
    model_name : str, optional (default="iJO1366")
        Model to test. Options: "iJO1366", "textbook", "e_coli_core"
    """
    import cobra
    import numpy as np
    from cobra.flux_analysis import moma, room
    from scipy.stats import pearsonr

    print("=" * 80)
    print(f"Testing Simulator vs COBRApy with {model_name} model")
    print("=" * 80)

    # Load model
    print(f"\n[1] Loading {model_name} model...")
    try:
        cobra_model = cobra.io.load_model(model_name)
        print(f"✓ Model loaded: {len(cobra_model.reactions)} reactions, {len(cobra_model.metabolites)} metabolites")
    except Exception as e:
        print(f"✗ Failed to load {model_name} model: {e}")
        print("  Trying alternative method...")
        try:
            cobra_model = cobra.io.read_sbml_model(f"{model_name}.xml")
            print("✓ Model loaded from file")
        except:
            print(f"✗ Please ensure {model_name} model is available")
            return

    # Initialize Simulator
    sim = Simulator()
    sim.load_cobra_model(cobra_model)

    # ===== FBA Test =====
    print("\n" + "=" * 80)
    print("[2] FBA (Flux Balance Analysis) Comparison")
    print("=" * 80)

    # Simulator FBA
    print("\nRunning Simulator FBA...")
    sim_status, sim_obj, sim_flux = sim.run_FBA()
    print(f"  Status: {sim_status}")
    print(f"  Objective value: {sim_obj:.6f}")

    # COBRApy FBA
    print("\nRunning COBRApy FBA...")
    cobra_solution = cobra_model.optimize()
    cobra_obj = cobra_solution.objective_value
    cobra_flux = cobra_solution.fluxes.to_dict()
    print(f"  Status: {cobra_solution.status}")
    print(f"  Objective value: {cobra_obj:.6f}")

    # Compare FBA results
    print("\nFBA Comparison:")
    obj_diff = abs(sim_obj - cobra_obj)
    print(f"  Objective difference: {obj_diff:.10f}")

    # Calculate correlation of flux distributions
    common_reactions = set(sim_flux.keys()) & set(cobra_flux.keys())
    sim_values = [sim_flux[r] for r in common_reactions]
    cobra_values = [cobra_flux[r] for r in common_reactions]

    if len(sim_values) > 1:
        corr, p_value = pearsonr(sim_values, cobra_values)
        print(f"  Flux correlation (Pearson): {corr:.10f}")

        # Calculate flux differences
        flux_diffs = [abs(sim_flux[r] - cobra_flux[r]) for r in common_reactions]
        max_diff = max(flux_diffs)
        mean_diff = np.mean(flux_diffs)
        print(f"  Max flux difference: {max_diff:.10f}")
        print(f"  Mean flux difference: {mean_diff:.10f}")

    # Store wild-type flux for MOMA/ROOM
    wild_flux = sim_flux.copy()

    # ===== Gene Knockout Setup =====
    print("\n" + "=" * 80)
    print("[3] Setting up gene knockout scenario")
    print("=" * 80)

    # Find a reaction to knockout (e.g., PGI - phosphoglucose isomerase)
    knockout_reaction = None
    for rxn_id in ["PGI", "PFK", "FBP"]:
        if rxn_id in sim.model_reactions:
            knockout_reaction = rxn_id
            break

    if knockout_reaction is None:
        # Just use the first reversible reaction
        knockout_reaction = sim.model_reactions[10]

    print(f"  Knocking out reaction: {knockout_reaction}")
    knockout_constraints = {knockout_reaction: (0, 0)}

    # Apply knockout to COBRApy model
    cobra_model_ko = cobra_model.copy()
    cobra_model_ko.reactions.get_by_id(knockout_reaction).bounds = (0, 0)

    # ===== Linear MOMA Test =====
    print("\n" + "=" * 80)
    print("[4] Linear MOMA (Minimization of Metabolic Adjustment) Comparison")
    print("=" * 80)

    # Simulator Linear MOMA
    print("\nRunning Simulator Linear MOMA...")
    sim_moma_status, sim_moma_obj, sim_moma_flux = sim.run_MOMA(
        wild_flux=wild_flux, flux_constraints=knockout_constraints
    )
    print(f"  Status: {sim_moma_status}")
    print(f"  Objective value (Manhattan distance): {sim_moma_obj:.6f}")

    # COBRApy Linear MOMA
    print("\nRunning COBRApy Linear MOMA...")
    try:
        cobra_moma_solution = moma(cobra_model_ko, solution=cobra_solution, linear=True)
        cobra_moma_obj = cobra_moma_solution.objective_value
        cobra_moma_flux = cobra_moma_solution.fluxes.to_dict()
        print(f"  Status: {cobra_moma_solution.status}")
        print(f"  Objective value: {cobra_moma_obj:.6f}")

        # Compare MOMA results
        print("\nLinear MOMA Comparison:")
        moma_obj_diff = abs(sim_moma_obj - cobra_moma_obj)
        print(f"  Objective difference: {moma_obj_diff:.10f}")

        # Flux correlation
        sim_moma_values = [sim_moma_flux[r] for r in common_reactions]
        cobra_moma_values = [cobra_moma_flux[r] for r in common_reactions]

        if len(sim_moma_values) > 1:
            moma_corr, _ = pearsonr(sim_moma_values, cobra_moma_values)
            print(f"  Flux correlation (Pearson): {moma_corr:.10f}")

            moma_flux_diffs = [abs(sim_moma_flux[r] - cobra_moma_flux[r]) for r in common_reactions]
            print(f"  Max flux difference: {max(moma_flux_diffs):.10f}")
            print(f"  Mean flux difference: {np.mean(moma_flux_diffs):.10f}")
    except Exception as e:
        print(f"  ✗ COBRApy MOMA failed: {e}")
        print("  Note: MOMA in COBRApy may require specific solvers")

    # ===== ROOM Test =====
    print("\n" + "=" * 80)
    print("[5] ROOM (Regulatory On/Off Minimization) Comparison")
    print("=" * 80)

    # Simulator ROOM
    print("\nRunning Simulator ROOM...")
    sim_room_status, sim_room_obj, sim_room_flux = sim.run_ROOM(
        wild_flux=wild_flux, flux_constraints=knockout_constraints, delta=0.03, epsilon=0.001
    )
    print(f"  Status: {sim_room_status}")
    print(f"  Number of changed reactions: {sim_room_obj:.0f}")

    # COBRApy ROOM
    print("\nRunning COBRApy ROOM...")
    try:
        cobra_room_solution = room(cobra_model_ko, solution=cobra_solution, delta=0.03, epsilon=0.001)
        cobra_room_obj = cobra_room_solution.objective_value
        cobra_room_flux = cobra_room_solution.fluxes.to_dict()
        print(f"  Status: {cobra_room_solution.status}")
        print(f"  Number of changed reactions: {cobra_room_obj:.0f}")

        # Compare ROOM results
        print("\nROOM Comparison:")
        room_obj_diff = abs(sim_room_obj - cobra_room_obj)
        print(f"  Changed reactions difference: {room_obj_diff:.0f}")

        # Flux correlation
        sim_room_values = [sim_room_flux[r] for r in common_reactions]
        cobra_room_values = [cobra_room_flux[r] for r in common_reactions]

        if len(sim_room_values) > 1:
            room_corr, _ = pearsonr(sim_room_values, cobra_room_values)
            print(f"  Flux correlation (Pearson): {room_corr:.10f}")

            room_flux_diffs = [abs(sim_room_flux[r] - cobra_room_flux[r]) for r in common_reactions]
            print(f"  Max flux difference: {max(room_flux_diffs):.10f}")
            print(f"  Mean flux difference: {np.mean(room_flux_diffs):.10f}")
    except Exception as e:
        print(f"  ✗ COBRApy ROOM failed: {e}")
        print("  Note: ROOM in COBRApy may require specific solvers")

    # ===== Summary =====
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"✓ FBA: Objective diff = {obj_diff:.10f}")
    print("✓ All tests completed")
    print("=" * 80)


if __name__ == "__main__":
    import sys

    # Parse command line arguments
    if len(sys.argv) > 1:
        test_type = sys.argv[1].lower()
        model = sys.argv[2] if len(sys.argv) > 2 else "textbook"
    else:
        test_type = "all"
        model = "textbook"

    if test_type in ["comparison", "all"]:
        print("\n" + "=" * 80)
        print("RUNNING SIMULATOR COMPARISON TESTS")
        print("=" * 80)
        test_simulator_comparison(model)

    # Usage instructions
    if len(sys.argv) == 1:
        print("\n" + "=" * 80)
        print("Usage: python simulator_linear_moma.py [test_type] [model]")
        print("  test_type: 'fseof', 'comparison', or 'all' (default: all)")
        print("  model: 'textbook', 'e_coli_core', 'iJO1366' (default: textbook)")
        print("=" * 80)
