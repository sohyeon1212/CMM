import Simulator
from gurobipy import *


class LAD(Simulator.Simulator):
    def __init__(self):
        """
        Constructor
        """

    def run_LP_fitting(self, opt_flux={}, flux_constraints={}, inf_flag=False):
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        Smatrix = self.Smatrix
        lower_bounds = self.lower_boundary_constraints.copy()
        upper_bounds = self.upper_boundary_constraints.copy()

        # Handle infinite bounds
        if not inf_flag:
            for key in lower_bounds:
                if lower_bounds[key] == float("-inf"):
                    lower_bounds[key] = -1000.0
            for key in upper_bounds:
                if upper_bounds[key] == float("inf"):
                    upper_bounds[key] = 1000.0

        pairs, coeffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model("RNASeq_DirectFlux")
        m.setParam("OutputFlag", 0)

        # Filter by opt_flux without applying a weight threshold
        target_reactions = {rid: w for rid, w in opt_flux.items() if rid in model_reactions and abs(w) > 0.01}

        if not target_reactions:
            print("No valid reactions above threshold.")
            return 1, None, None

        # Variable creation
        v, f, b, delta = {}, {}, {}, {}
        for reaction_id in model_reactions:
            # Apply flux constraints first
            if reaction_id in flux_constraints:
                lb, ub = flux_constraints[reaction_id]
            else:
                lb, ub = lower_bounds[reaction_id], upper_bounds[reaction_id]

            v[reaction_id] = m.addVar(lb=lb, ub=ub, name=f"v_{reaction_id}")
            f[reaction_id] = m.addVar(lb=0.0, name=f"f_{reaction_id}")
            b[reaction_id] = m.addVar(lb=0.0, name=f"b_{reaction_id}")

        for reaction_id in target_reactions:
            delta[reaction_id] = m.addVar(lb=0.0, name=f"delta_{reaction_id}")

        m.update()

        # Flux decomposition and irreversibility constraints
        for reaction_id in model_reactions:
            m.addConstr(v[reaction_id] == f[reaction_id] - b[reaction_id])
            if lower_bounds[reaction_id] >= 0:
                m.addConstr(b[reaction_id] == 0)
            if upper_bounds[reaction_id] <= 0:
                m.addConstr(f[reaction_id] == 0)

        # Stoichiometric constraints
        for met_id in model_metabolites:
            if len(pairs.select(met_id, "*")) > 0:
                m.addConstr(
                    quicksum(v[rxn_id] * coeffvalue[met_id, rxn_id] for met_id, rxn_id in pairs.select(met_id, "*"))
                    == 0
                )

        # Absolute deviation constraints
        scaling_factor = 1.0
        epsilon = 1e-9
        for reaction_id, weight in target_reactions.items():
            target_flux = abs(weight) * scaling_factor
            m.addConstr(delta[reaction_id] >= (f[reaction_id] + b[reaction_id]) - target_flux)
            m.addConstr(delta[reaction_id] >= target_flux - (f[reaction_id] + b[reaction_id]))

        m.update()

        objective = quicksum(delta[rid] for rid in target_reactions)
        m.setObjective(objective, GRB.MINIMIZE)

        m.optimize()

        if m.status == GRB.Status.OPTIMAL:
            flux_results = {rid: v[rid].x for rid in model_reactions}
            return m.status, m.ObjVal, flux_results
        else:
            return m.status, None, None


def read_expression_data(filename):
    expression_info = {}
    fp = open(filename)
    fp.readline()
    for line in fp:
        sptlist = line.split("\t")
        gene_id = sptlist[0].strip()
        value = sptlist[1].strip()
        expression_info[gene_id] = float(value)
    fp.close()
    return expression_info
