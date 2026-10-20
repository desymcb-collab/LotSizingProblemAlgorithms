#include "WSN_MODEL.h"

using namespace std;

WSN_MODEL::WSN_MODEL(const IloEnv env, const WSNInput* input, const char* modelName) : env(env), model(env),
			 adj(*(input->adj)), e0(*(input->e0)), d(*(input->d)),
			 x(env), b(env), e(env), r(env), nOut(env), nIn(env)
{
//  strcpy(this->modelName,modelName);

  try
  {
    this->buildInput(input);

    this->buildVariables();

    this->buildConstraints();
  }

  catch(std::bad_alloc& ba)
  {
    env.out() << "Bad allocation of " << modelName << " \n";

    throw ba;
  }
}

// ========================================================================

void WSN_MODEL::buildInput(const WSNInput* input)
{
  	this->N = input->N;
  	this->T = input->T;
}

void WSN_MODEL::buildVariables()
{
  this->env.out() << "Initializing flow variables and expressions\n";

  // Initialize a flow array for every time interval
  this->x.add(T, IloBoolVarArray(this->env));

  // Initialize a neighborhood for every time interval
  this->nOut.add(T, IloArray<IloExpr> (this->env));
  this->nIn.add(T, IloArray<IloExpr> (this->env));

  // Iterate through time intervals
  for (IloInt t = 0; t < this->T; t++)
  {
    // Initialize a flow array for every t
    this->x[t] = IloBoolVarArray(this->env);

    // Initialize a neighborhood for every t
    this->nOut[t] = IloArray<IloExpr> (this->env, this->N);
    this->nIn[t] = IloArray<IloExpr> (this->env, this->N);

    // Initialize a neighborhood expression for every node
    for (IloInt i = 0; i < this->N; i++)
    {
      this->nOut[t][i] = IloExpr(this->env);
      this->nIn[t][i] = IloExpr(this->env);
  	}

    IloInt arc = 0;

    // Iterate through arcs
    for (IloInt i = 0; i < this->N; i++)
    {
      for (IloInt j = 0; j < this->N; j++)
      {
        // Check arc feasibility
        if (this->adj[i][j])
        {
          //sprintf(xName, "x({%ld,%ld})(%ld)", i, j, t); x[t].add(IloBoolVar(env, xName));

          // Add a flow variable for every time inerval and feasible arc
          this->x[t].add(IloBoolVar(this->env));

          // Update its neighborhood
          this->nOut[t][i] += this->x[t][arc];
          this->nIn[t][j] += this->x[t][arc++];
        }
      }
    }
  }

  // ==========================================================================

  this->env.out() << "Initializing energy, buffer and rx variables\n";

  // Initialize an vararray for every node
  b.add(this->N, IloIntVarArray(this->env));
  e.add(this->N, IloIntVarArray(this->env));
  r.add(this->N, IloIntVarArray(this->env));

  // Iterate through nodes
  for (IloInt i = 0; i < this->N; i++)
  {
    // Initialize an vararray for every node
    b[i] = IloIntVarArray(this->env);
    e[i] = IloIntVarArray(this->env);
    r[i] = IloIntVarArray(this->env);

    // Iterate through time intervals
    for (IloInt t = 0; t < this->T; t++)
    {
      //			sprintf(xName, "b(%ld)(%ld)", i, t); b[i].add(IloIntVar(env, xName));
      //			sprintf(xName, "_e(%ld)(%ld)", i, t); e[i].add(IloIntVar(env, xName));
      //			sprintf(xName, "Rx(%ld)(%ld)", i, t); r[i].add(IloIntVar(env, xName));

      // Add a variable for every node and time interval
      b[i].add(IloIntVar(this->env));
      e[i].add(IloIntVar(this->env));
      r[i].add(IloIntVar(this->env));
    }
  }
}

// ========================================================================

void WSN_MODEL::buildEnergyConstraints()
{
  env.out() << "Building Energy Constraints...\n";

  // Iterate over non-sink nodes
	for (IloInt i = 1; i < N; i++)
	{
    // Initial Constraint
    model.add(e[i][0] == e0[i] - nOut[0][i]);

    // Iterate over time intervals
		for (IloInt t = 1; t < T; t++)
			model.add(e[i][t] == e[i][t-1] - nOut[t][i]);
	}
}

void WSN_MODEL::buildSupplyDemandConstraints()
{
  env.out() << "Building Supply/Demand Constraints...\n";

  // Iterate over nodes
	for (IloInt i = 0; i < N; i++)
	{
    // Initial Constraint
    model.add(b[i][0] == d[i][0] - nOut[0][i]);

    // Iterate over time intervals
		for (IloInt t = 1; t < T; t++)
			model.add(b[i][t] == b[i][t-1] + d[i][t] + r[i][t] - nOut[t][i]);
	}
}

void WSN_MODEL::buildTxConstraints()
{
  env.out() << "Building Sending Constraints...\n";

  // Iterate over edges not originated from sink-nodes
  for (IloInt i = 1; i < N; i++)
  {
    for (IloInt j = 0; j < N; j++)
    {
      // Check if edge is feasible
      if (adj[i][j])
      {
        // Iterate over time intervals
        for (IloInt t = 0; t < T; t++)
        {
          model.add(nOut[t][i] + nOut[t][j] <= 1);
          // model.add(IloIfThen(env, r[j][t] + d[j][t] > 0, nOut[t][i] + nOut[t][j] > 0));
        }
      }

    }
  }
}

void WSN_MODEL::buildRxConstraints()
{
  env.out() << "Building Reception Constraints...\n";

  // Iterate over time intervals
  for (IloInt t = 1; t < T; t++)
  {
    // Iterate over nodes
    for (IloInt i = 0; i < N; i++)
    {
      model.add(IloIfThen(env, nIn[t-1][i] > 1, r[i][t] == 0));

      // model.add(r[i][t] <= nIn[t-1][i]);
			model.add(IloIfThen(env, nIn[t-1][i] <= 1, r[i][t] == nIn[t-1][i]));
    }
  }
}

void WSN_MODEL::buildConstraints()
{
  this->buildEnergyConstraints();
  this->buildSupplyDemandConstraints();
  this->buildTxConstraints();
  this->buildRxConstraints();
}

// ========================================================================

IloObjective WSN_MODEL::buildEnergyObjective()
{
  IloExpr energy(env);

  // Iterate over non-sink nodes
	for (IloInt i = 1; i < N; i++)
		energy += e[i][T-1];

//	IloExtractable energy_obj = model.add(IloMaximize(env, energy));

//	return energy_obj;

	return IloMaximize(env, energy);
}

IloObjective WSN_MODEL::buildDemandObjective()
{
//  IloExtractable demand_obj = model.add(IloMaximize(env, b[0][T-1]));

//	return demand_obj;

	return IloMaximize(env, b[0][T-1]);
}
