// generated with the command : ../tools/builder -I/usr/local/include/Siconos/Kernel -I/usr/local/include/Siconos/Numerics -I/usr/include/libxml2
SICONOS_IO_REGISTER_WITH_BASES(ControlFirstOrderLinearDS, (ControlDynamicalSystem),
                               (_x0)
                               (_A))
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerFrom1DLocalFrameR, (NewtonEulerR),
                               (_isOnContact)
                               (_Pc1)
                               (_Pc2)
                               (_Nc)
                               (_Mabs_C)
                               (_NPG1)
                               (_NPG2)
                               (_AUX1)
                               (_AUX2))
SICONOS_IO_REGISTER_WITH_BASES(NonSmoothEvent, (Event),
                              )
SICONOS_IO_REGISTER_WITH_BASES(OSNSMatrixProjectOnConstraints, (OSNSMatrix),
                              )
SICONOS_IO_REGISTER_WITH_BASES(BlockMatrix, (SiconosMatrix),
                               (_mat)
                               (_tabRow)
                               (_tabCol))
SICONOS_IO_REGISTER_WITH_BASES(BlockVector, (SiconosVector),
                               (_sizeV)
                               (vect)
                               (_tabIndex))
SICONOS_IO_REGISTER(NonSmoothDynamicalSystem,
                    (_BVP)
                    (_topology)
                    (_mIsLinear))
SICONOS_IO_REGISTER(Relation,
                    (_pluginh)
                    (_pluginJachx)
                    (_pluginJachlambda)
                    (_pluging)
                    (_pluginJacLg)
                    (_pluginf)
                    (_plugine)
                    (relationType)
                    (subType)
                    (_interaction)
                    (data)
                    (_workR)
                    (_workX)
                    (_workXdot)
                    (_workZ)
                    (_workY)
                    (_workL)
                    (_Residuy)
                    (_h_alpha)
                    (_jachlambda))
SICONOS_IO_REGISTER(ControlManager,
                    (_allSensors)
                    (_allActuators)
                    (_model))
SICONOS_IO_REGISTER_WITH_BASES(NewtonImpactNSL, (NonSmoothLaw),
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(LinearSMC, (CommonSMC),
                               (_SMC)
                               (_DS_SMC)
                               (_tD_SMC)
                               (_simulationSMC)
                               (_integratorSMC)
                               (_thetaSMC)
                               (_OSNSPB_SMC)
                               (_sampledControl)
                               (_eventsManager)
                               (_nsLawSMC))
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerDS, (DynamicalSystem),
                               (_v)
                               (_v0)
                               (_vMemory)
                               (_qMemory)
                               (_forcesMemory)
                               (_dotqMemory)
                               (_qDim)
                               (_q)
                               (_deltaq)
                               (_q0)
                               (_dotq)
                               (_MObjToAbs)
                               (_I)
                               (_mass)
                               (_massMatrix)
                               (_luW)
                               (_T)
                               (_p)
                               (_mExt)
                               (_fExt)
                               (_forces)
                               (_jacobianvFL)
                               (_jacobianqDotForces))
SICONOS_IO_REGISTER_WITH_BASES(ControlSensor, (Sensor),
                               (_YDim)
                               (_storedY))
SICONOS_IO_REGISTER(OSNSMatrix,
                    (dimRow)
                    (dimColumn)
                    (storageType)
                    (DSBlocksPositions)
                    (M1)
                    (Mt)
                    (M2))
SICONOS_IO_REGISTER_WITH_BASES(RelayNSL, (NonSmoothLaw),
                               (_lb)
                               (_ub))
SICONOS_IO_REGISTER_WITH_BASES(MixedComplementarityConditionNSL, (NonSmoothLaw),
                               (EqualitySize))
SICONOS_IO_REGISTER_WITH_BASES(CommonSMC, (Actuator),
                               (_sDim)
                               (_indx)
                               (_u)
                               (_Csurface)
                               (_sensor)
                               (_initDone)
                               (_curDeltaT)
                               (_B)
                               (_D)
                               (_relationSMC)
                               (_sign)
                               (_interactionSMC)
                               (_lambda)
                               (_xController)
                               (_precision))
SICONOS_IO_REGISTER(TimeDiscretisation,
                    (_h)
                    (_k)
                    (_tk)
                    (_tdCase)
                    (_pos))
SICONOS_IO_REGISTER_WITH_BASES(SensorEvent, (Event),
                               (_sensor))
SICONOS_IO_REGISTER(OneStepIntegrator,
                    (integratorType)
                    (OSIDynamicalSystems)
                    (OSIInteractions)
                    (_sizeMem)
                    (simulationLink))
SICONOS_IO_REGISTER_WITH_BASES(MLCP, (LinearOSNS),
                               (_n)
                               (_m)
                               (_curBlock))
SICONOS_IO_REGISTER_WITH_BASES(TimeSteppingD1Minus, (Simulation),
                               (impactOccuredLastTimeStep))
SICONOS_IO_REGISTER(GraphProperties,
                    (symmetric))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianR, (Relation),
                               (_jachq)
                               (_jachqDot))
SICONOS_IO_REGISTER_WITH_BASES(SchatzmanPaoli, (OneStepIntegrator),
                               (WMap)
                               (_WBoundaryConditionsMap)
                               (_theta)
                               (_gamma)
                               (_useGamma)
                               (_useGammaForRelation))
SICONOS_IO_REGISTER_WITH_BASES(LinearChatteringSMC, (CommonSMC),
                               (_s))
SICONOS_IO_REGISTER(SiconosMatrix,
                    (dimRow)
                    (dimCol)
                    (num))
SICONOS_IO_REGISTER(BlockCSRMatrix,
                    (nr)
                    (nc)
                    (diagSizes)
                    (rowPos)
                    (colPos))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianLinearTIR, (LagrangianR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(Moreau2, (Moreau),
                              )
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerFrom3DLocalFrameR, (NewtonEulerFrom1DLocalFrameR),
                              )
SICONOS_IO_REGISTER_WITH_BASES(LinearSensor, (ControlSensor),
                               (_data)
                               (_dataPlot)
                               (_k)
                               (_matC)
                               (_matD)
                               (_nSteps))
SICONOS_IO_REGISTER_WITH_BASES(SampledPIDActuator, (Actuator),
                               (_ref)
                               (_u)
                               (_K)
                               (_sensor)
                               (_initDone)
                               (_curDeltaT))
SICONOS_IO_REGISTER_WITH_BASES(LCP, (LinearOSNS),
                              )
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerR, (Relation),
                               (_ysize)
                               (_xsize)
                               (_qsize)
                               (_workQ)
                               (_jachq)
                               (_jachqDot)
                               (_jachlambda)
                               (_e)
                               (_contactForce)
                               (_jachqT))
SICONOS_IO_REGISTER_WITH_BASES(EventDriven, (Simulation),
                               (_istate)
                               (TOL_ED))
SICONOS_IO_REGISTER(ControlDynamicalSystem,
                    (_t0)
                    (_T)
                    (_h)
                    (_theta)
                    (_elapsedTime)
                    (_N)
                    (_nDim)
                    (_x0)
                    (_dataM)
                    (_processDS)
                    (_model)
                    (_processTD)
                    (_processSimulation)
                    (_processIntegrator)
                    (_CM))
SICONOS_IO_REGISTER(OneStepNSProblem,
                    (_id)
                    (_sizeOutput)
                    (_DSBlocks)
                    (_simulation)
                    (_OSNSInteractions)
                    (_levelMin)
                    (_levelMax)
                    (_maxSize)
                    (_CPUtime)
                    (_nbIter)
                    (_hasBeenUpdated))
SICONOS_IO_REGISTER_WITH_BASES(TimeStepping, (Simulation),
                               (_computeResiduY)
                               (_computeResiduR)
                               (_newtonTolerance)
                               (_newtonMaxIteration)
                               (_newtonNbSteps)
                               (_newtonOptions)
                               (_newtonResiduDSMax)
                               (_newtonResiduYMax)
                               (_newtonResiduRMax))
SICONOS_IO_REGISTER(Interaction,
                    (_initialized)
                    (_id)
                    (_number)
                    (_relativeDegree)
                    (_lowerLevelForOutput)
                    (_upperLevelForOutput)
                    (_lowerLevelForInput)
                    (_upperLevelForInput)
                    (_interactionSize)
                    (_sizeOfDS)
                    (_sizeZ)
                    (_absolutePosition)
                    (_absolutePositionProj)
                    (_y)
                    (_yOld)
                    (_y_k)
                    (_yMemory)
                    (_steps)
                    (_lambda)
                    (_lambdaOld)
                    (_involvedDS)
                    (_nslaw)
                    (_relation)
                    (_workX)
                    (_workXq)
                    (_workFree)
                    (_workYp)
                    (_workZ))
SICONOS_IO_REGISTER(DynamicalSystem,
                    (_number)
                    (_n)
                    (_x0)
                    (_residuFree)
                    (_r)
                    (_normRef)
                    (_x)
                    (_jacxRhs)
                    (_jacgx)
                    (_jacxDotG)
                    (_z)
                    (_g)
                    (_pluging)
                    (_pluginJacgx)
                    (_pluginJacxDotG)
                    (_xMemory)
                    (_stepsInMemory)
                    (_workV)
                    (_workMatrix)
                    (_workFree)
                    (count))
SICONOS_IO_REGISTER(Sensor,
                    (_type)
                    (_id)
                    (_model)
                    (_DS)
                    (_DSx)
                    (_nDim)
                    (_timeDiscretisation)
                    (_eSensor))
SICONOS_IO_REGISTER(Topology,
                    (_allInteractions)
                    (_DSG)
                    (_IG)
                    (_isTopologyUpToDate)
                    (_hasChanged)
                    (_numberOfConstraints)
                    (_symmetric))
SICONOS_IO_REGISTER(Actuator,
                    (_type)
                    (_nDim)
                    (_id)
                    (_allSensors)
                    (_allDS)
                    (_DS)
                    (_model)
                    (_timeDiscretisation)
                    (_eActuator))
SICONOS_IO_REGISTER(SiconosException,
                    (reportMsg))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearTIDS, (FirstOrderLinearDS),
                              )
SICONOS_IO_REGISTER(SiconosMemory,
                    (_maxSize)
                    (_nbVectorsInMemory)
                    (_vectorMemory))
SICONOS_IO_REGISTER(InteractionProperties,
                    (block)
                    (source)
                    (target)
                    (forControl))
SICONOS_IO_REGISTER(Simulation,
                    (_name)
                    (_timeDiscretisation)
                    (_eventsManager)
                    (_tinit)
                    (_tend)
                    (_tout)
                    (_allOSI)
                    (_osiMap)
                    (_interactionOsiMap)
                    (_allNSProblems)
                    (_model)
                    (_levelMinForOutput)
                    (_levelMaxForOutput)
                    (_levelMinForInput)
                    (_levelMaxForInput)
                    (_tolerance)
                    (_printStat)
                    (_staticLevels)
                    (_levelsAreComputed)
                    (statOut)
                    (_useRelativeConvergenceCriterion)
                    (_relativeConvergenceCriterionHeld)
                    (_relativeConvergenceTol))
SICONOS_IO_REGISTER(SystemProperties,
                    (upper_block)
                    (lower_block))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianLinearTIDS, (LagrangianDS),
                               (_K)
                               (_C))
SICONOS_IO_REGISTER_WITH_BASES(GenericMechanical, (LinearOSNS),
                              )
SICONOS_IO_REGISTER_WITH_BASES(LagrangianScleronomousR, (LagrangianR),
                               (_pluginjqh)
                               (_pluginjqhdot)
                               (_NLh2dot))
SICONOS_IO_REGISTER(Model,
                    (_t)
                    (_t0)
                    (_T)
                    (_strat)
                    (_nsds)
                    (_title)
                    (_author)
                    (_description)
                    (_date))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderNonLinearDS, (DynamicalSystem),
                               (_M)
                               (_f)
                               (_fold)
                               (_jacobianfx)
                               (_pluginf)
                               (_pluginJacxf)
                               (_pluginM)
                               (_rMemory)
                               (_residur)
                               (_g_alpha)
                               (_xp)
                               (_xq)
                               (_invM))
SICONOS_IO_REGISTER_WITH_BASES(Relay, (LinearOSNS),
                               (_lb)
                               (_ub))
SICONOS_IO_REGISTER_WITH_BASES(ComplementarityConditionNSL, (NonSmoothLaw),
                              )
SICONOS_IO_REGISTER_WITH_BASES(EqualityConditionNSL, (NonSmoothLaw),
                              )
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearDS, (FirstOrderNonLinearDS),
                               (_A)
                               (_b)
                               (_pluginA)
                               (_pluginb))
SICONOS_IO_REGISTER_WITH_BASES(LinearOSNS, (OneStepNSProblem),
                               (_w)
                               (_z)
                               (_M)
                               (_q)
                               (_MStorageType)
                               (_keepLambdaAndYState))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderType2R, (FirstOrderR),
                               (jacgx))
SICONOS_IO_REGISTER_WITH_BASES(TimeSteppingProjectOnConstraints, (TimeStepping),
                               (_indexSetLevelForProjection)
                               (_constraintTol)
                               (_constraintTolUnilateral)
                               (_projectionMaxIteration)
                               (_doProj)
                               (_doOnlyProj))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianRheonomousR, (LagrangianR),
                               (_hDot)
                               (_pluginhDot)
                               (_pluginJachq))
SICONOS_IO_REGISTER_WITH_BASES(TimeDiscretisationEvent, (Event),
                              )
SICONOS_IO_REGISTER_WITH_BASES(MultipleImpactNSL, (NonSmoothLaw),
                               (_ResCof)
                               (_Stiff)
                               (_ElasCof))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianCompliantR, (LagrangianR),
                               (_pluginJachq)
                               (_pluginJachlambda))
SICONOS_IO_REGISTER_WITH_BASES(TimeSteppingCombinedProjection, (TimeStepping),
                               (_indexSetLevelForProjection)
                               (_constraintTol)
                               (_constraintTolUnilateral)
                               (_projectionMaxIteration)
                               (_doCombinedProj)
                               (_isIndexSetsStable))
SICONOS_IO_REGISTER(NonSmoothLaw,
                    (_size)
                    (_sizeProjectOnConstraints))
SICONOS_IO_REGISTER(BoundaryCondition,
                    (_velocityIndices)
                    (_prescribedVelocity)
                    (_prescribedVelocityOld)
                    (_pluginPrescribedVelocity))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearR, (FirstOrderR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(InteractionsGraph, (_InteractionsGraph),
                               (blockProj)
                               (upper_blockProj)
                               (lower_blockProj)
                               (dummy))
SICONOS_IO_REGISTER_WITH_BASES(MoreauCombinedProjectionOSI, (Moreau),
                              )
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearTIR, (FirstOrderR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(Equality, (LinearOSNS),
                              )
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderR, (Relation),
                               (Jachx)
                               (Jacglambda))
SICONOS_IO_REGISTER_WITH_BASES(Moreau, (OneStepIntegrator),
                               (WMap)
                               (_WBoundaryConditionsMap)
                               (_theta)
                               (_gamma)
                               (_useGamma)
                               (_useGammaForRelation))
SICONOS_IO_REGISTER(Event,
                    (timeOfEvent)
                    (type)
                    (dTime)
                    (tick))
SICONOS_IO_REGISTER_WITH_BASES(ActuatorEvent, (Event),
                               (_actuator))
SICONOS_IO_REGISTER_WITH_BASES(DynamicalSystemsGraph, (_DynamicalSystemsGraph),
                               (OSI)
                               (dummy))
SICONOS_IO_REGISTER_WITH_BASES(MLCPProjectOnConstraints, (MLCP),
                               (_alpha))
SICONOS_IO_REGISTER_WITH_BASES(OSNSMultipleImpact, (LinearOSNS),
                               (Impulse_variable)
                               (Time_variable)
                               (Ncontact)
                               (NstepEst)
                               (NstepMax)
                               (TOL_IMPACT)
                               (TypeCompLaw)
                               (VelContact)
                               (OldVelContact)
                               (EnerContact)
                               (WcContact)
                               (DistriVector)
                               (StateContact)
                               (Kcontact)
                               (ResContact)
                               (ElasCoefContact)
                               (DelImpulseContact)
                               (TolImpulseContact)
                               (ImpulseContact_update)
                               (ForceContact)
                               (SelectPrimaConInVel)
                               (IdPrimaContact)
                               (IsPrimaConEnergy)
                               (VelAtPrimaCon)
                               (EnerAtPrimaCon)
                               (DeltaP)
                               (OutputFile)
                               (NameFile)
                               (YesSaveData)
                               (NstepSave)
                               (IsNumberOfStepsEst)
                               (_DataMatrix)
                               (YesSaveByMatrix)
                               (SizeDataSave)
                               (_IsImpactEnd))
SICONOS_IO_REGISTER(PluggedObject,
                    (_pluginName))
SICONOS_IO_REGISTER_WITH_BASES(NewtonImpactFrictionNSL, (NonSmoothLaw),
                               (_en)
                               (_et)
                               (_mu))
SICONOS_IO_REGISTER_WITH_BASES(QP, (OneStepNSProblem),
                               (Q)
                               (_p))
SICONOS_IO_REGISTER(EventsManager,
                    (_allEvents)
                    (_currentEvent)
                    (_nextEvent)
                    (_ETD)
                    (_ENonSmooth)
                    (_simulation)
                    (_hasNS)
                    (_hasCM)
                    (_GapLimit2Events))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianDS, (DynamicalSystem),
                               (_ndof)
                               (_q)
                               (_q0)
                               (_velocity0)
                               (_qMemory)
                               (_velocityMemory)
                               (_p)
                               (_mass)
                               (_fInt)
                               (_jacobianFIntq)
                               (_jacobianFIntqDot)
                               (_fExt)
                               (_NNL)
                               (_jacobianNNLq)
                               (_jacobianNNLqDot)
                               (_forces)
                               (_jacobianqForces)
                               (_jacobianqDotForces)
                               (_boundaryConditions)
                               (_reactionToBoundaryConditions)
                               (_pluginMass)
                               (_pluginFInt)
                               (_pluginFExt)
                               (_pluginNNL)
                               (_pluginJacqFInt)
                               (_pluginJacqDotFInt)
                               (_pluginJacqNNL)
                               (_pluginJacqDotNNL))
SICONOS_IO_REGISTER_WITH_BASES(Circle, (CircularDS),
                              )
SICONOS_IO_REGISTER_WITH_BASES(CircularDS, (LagrangianDS),
                               (radius)
                               (massValue))
SICONOS_IO_REGISTER_WITH_BASES(CircleCircleR, (CircularR),
                               (ar1mr2))
SICONOS_IO_REGISTER_WITH_BASES(CircularR, (LagrangianScleronomousR),
                               (r1)
                               (r2))
SICONOS_IO_REGISTER_WITH_BASES(Disk, (CircularDS),
                              )
SICONOS_IO_REGISTER_WITH_BASES(DiskDiskR, (CircularR),
                               (r1pr2))
SICONOS_IO_REGISTER_WITH_BASES(DiskMovingPlanR, (LagrangianRheonomousR),
                               (_time)
                               (_A)
                               (_B)
                               (_C)
                               (_ADot)
                               (_BDot)
                               (_CDot)
                               (_sqrA2pB2)
                               (_r)
                               (_AADot)
                               (_BBDot)
                               (_cubsqrA2pB2)
                               (_AFunction)
                               (_BFunction)
                               (_CFunction)
                               (_ADotFunction)
                               (_BDotFunction)
                               (_CDotFunction))
SICONOS_IO_REGISTER_WITH_BASES(DiskPlanR, (LagrangianScleronomousR),
                               (r)
                               (A)
                               (B)
                               (C)
                               (sqrA2pB2)
                               (AC)
                               (B2)
                               (A2)
                               (AB)
                               (BC)
                               (xCenter)
                               (yCenter)
                               (width)
                               (halfWidth)
                               (x1)
                               (x2)
                               (y1)
                               (y2)
                               (finite))
SICONOS_IO_REGISTER_WITH_BASES(SphereLDS, (LagrangianDS),
                               (radius)
                               (massValue)
                               (I))
SICONOS_IO_REGISTER_WITH_BASES(SphereLDSPlanR, (LagrangianScleronomousR),
                               (r)
                               (A)
                               (B)
                               (C)
                               (D)
                               (nN)
                               (nU)
                               (u1)
                               (u2)
                               (u3)
                               (v1)
                               (v2)
                               (v3)
                               (n1)
                               (n2)
                               (n3)
                               (ru1)
                               (ru2)
                               (ru3)
                               (rv1)
                               (rv2)
                               (rv3))
SICONOS_IO_REGISTER_WITH_BASES(SphereLDSSphereLDSR, (LagrangianScleronomousR),
                               (r1)
                               (r2)
                               (r1pr2))
SICONOS_IO_REGISTER_WITH_BASES(SphereNEDS, (NewtonEulerDS),
                               (radius))
SICONOS_IO_REGISTER_WITH_BASES(SphereNEDSPlanR, (NewtonEulerFrom3DLocalFrameR),
                               (r)
                               (A)
                               (B)
                               (C)
                               (D)
                               (nN)
                               (nU)
                               (u1)
                               (u2)
                               (u3)
                               (v1)
                               (v2)
                               (v3)
                               (n1)
                               (n2)
                               (n3)
                               (ru1)
                               (ru2)
                               (ru3)
                               (rv1)
                               (rv2)
                               (rv3))
SICONOS_IO_REGISTER_WITH_BASES(SphereNEDSSphereNEDSR, (NewtonEulerFrom3DLocalFrameR),
                               (r1)
                               (r2)
                               (r1pr2))
SICONOS_IO_REGISTER(SiconosBodies,
                    (_plans)
                    (_model)
                    (_playground))
SICONOS_IO_REGISTER(Hashed,
                    (body)
                    (i)
                    (j)
                    (k))
SICONOS_IO_REGISTER(SpaceFilter,
                    (_bboxfactor)
                    (_cellsize)
                    (_interID)
                    (_model)
                    (_nslaw)
                    (_plans)
                    (_hash_table)
                    (_osnsinit))

template <class Archive>
void siconos_io_register_generated(Archive& ar)
{
  ar.register_type(static_cast<ControlFirstOrderLinearDS*>(NULL));
  ar.register_type(static_cast<NewtonEulerFrom1DLocalFrameR*>(NULL));
  ar.register_type(static_cast<NonSmoothEvent*>(NULL));
  ar.register_type(static_cast<OSNSMatrixProjectOnConstraints*>(NULL));
  ar.register_type(static_cast<BlockMatrix*>(NULL));
  ar.register_type(static_cast<BlockVector*>(NULL));
  ar.register_type(static_cast<NonSmoothDynamicalSystem*>(NULL));
  ar.register_type(static_cast<ControlManager*>(NULL));
  ar.register_type(static_cast<NewtonImpactNSL*>(NULL));
  ar.register_type(static_cast<LinearSMC*>(NULL));
  ar.register_type(static_cast<NewtonEulerDS*>(NULL));
  ar.register_type(static_cast<OSNSMatrix*>(NULL));
  ar.register_type(static_cast<RelayNSL*>(NULL));
  ar.register_type(static_cast<MixedComplementarityConditionNSL*>(NULL));
  ar.register_type(static_cast<TimeDiscretisation*>(NULL));
  ar.register_type(static_cast<SensorEvent*>(NULL));
  ar.register_type(static_cast<MLCP*>(NULL));
  ar.register_type(static_cast<TimeSteppingD1Minus*>(NULL));
  ar.register_type(static_cast<GraphProperties*>(NULL));
  ar.register_type(static_cast<SchatzmanPaoli*>(NULL));
  ar.register_type(static_cast<LinearChatteringSMC*>(NULL));
  ar.register_type(static_cast<BlockCSRMatrix*>(NULL));
  ar.register_type(static_cast<LagrangianLinearTIR*>(NULL));
  ar.register_type(static_cast<Moreau2*>(NULL));
  ar.register_type(static_cast<NewtonEulerFrom3DLocalFrameR*>(NULL));
  ar.register_type(static_cast<LinearSensor*>(NULL));
  ar.register_type(static_cast<SampledPIDActuator*>(NULL));
  ar.register_type(static_cast<LCP*>(NULL));
  ar.register_type(static_cast<NewtonEulerR*>(NULL));
  ar.register_type(static_cast<EventDriven*>(NULL));
  ar.register_type(static_cast<ControlDynamicalSystem*>(NULL));
  ar.register_type(static_cast<TimeStepping*>(NULL));
  ar.register_type(static_cast<Interaction*>(NULL));
  ar.register_type(static_cast<Topology*>(NULL));
  ar.register_type(static_cast<SiconosException*>(NULL));
  ar.register_type(static_cast<FirstOrderLinearTIDS*>(NULL));
  ar.register_type(static_cast<SiconosMemory*>(NULL));
  ar.register_type(static_cast<InteractionProperties*>(NULL));
  ar.register_type(static_cast<SystemProperties*>(NULL));
  ar.register_type(static_cast<LagrangianLinearTIDS*>(NULL));
  ar.register_type(static_cast<GenericMechanical*>(NULL));
  ar.register_type(static_cast<LagrangianScleronomousR*>(NULL));
  ar.register_type(static_cast<Model*>(NULL));
  ar.register_type(static_cast<FirstOrderNonLinearDS*>(NULL));
  ar.register_type(static_cast<Relay*>(NULL));
  ar.register_type(static_cast<ComplementarityConditionNSL*>(NULL));
  ar.register_type(static_cast<EqualityConditionNSL*>(NULL));
  ar.register_type(static_cast<FirstOrderLinearDS*>(NULL));
  ar.register_type(static_cast<FirstOrderType2R*>(NULL));
  ar.register_type(static_cast<TimeSteppingProjectOnConstraints*>(NULL));
  ar.register_type(static_cast<LagrangianRheonomousR*>(NULL));
  ar.register_type(static_cast<TimeDiscretisationEvent*>(NULL));
  ar.register_type(static_cast<MultipleImpactNSL*>(NULL));
  ar.register_type(static_cast<LagrangianCompliantR*>(NULL));
  ar.register_type(static_cast<TimeSteppingCombinedProjection*>(NULL));
  ar.register_type(static_cast<BoundaryCondition*>(NULL));
  ar.register_type(static_cast<FirstOrderLinearR*>(NULL));
  ar.register_type(static_cast<InteractionsGraph*>(NULL));
  ar.register_type(static_cast<MoreauCombinedProjectionOSI*>(NULL));
  ar.register_type(static_cast<FirstOrderLinearTIR*>(NULL));
  ar.register_type(static_cast<Equality*>(NULL));
  ar.register_type(static_cast<Moreau*>(NULL));
  ar.register_type(static_cast<ActuatorEvent*>(NULL));
  ar.register_type(static_cast<DynamicalSystemsGraph*>(NULL));
  ar.register_type(static_cast<MLCPProjectOnConstraints*>(NULL));
  ar.register_type(static_cast<OSNSMultipleImpact*>(NULL));
  ar.register_type(static_cast<PluggedObject*>(NULL));
  ar.register_type(static_cast<NewtonImpactFrictionNSL*>(NULL));
  ar.register_type(static_cast<EventsManager*>(NULL));
  ar.register_type(static_cast<LagrangianDS*>(NULL));
  ar.register_type(static_cast<Circle*>(NULL));
  ar.register_type(static_cast<CircularDS*>(NULL));
  ar.register_type(static_cast<CircleCircleR*>(NULL));
  ar.register_type(static_cast<CircularR*>(NULL));
  ar.register_type(static_cast<Disk*>(NULL));
  ar.register_type(static_cast<DiskDiskR*>(NULL));
  ar.register_type(static_cast<DiskMovingPlanR*>(NULL));
  ar.register_type(static_cast<DiskPlanR*>(NULL));
  ar.register_type(static_cast<SphereLDS*>(NULL));
  ar.register_type(static_cast<SphereLDSPlanR*>(NULL));
  ar.register_type(static_cast<SphereLDSSphereLDSR*>(NULL));
  ar.register_type(static_cast<SphereNEDS*>(NULL));
  ar.register_type(static_cast<SphereNEDSPlanR*>(NULL));
  ar.register_type(static_cast<SphereNEDSSphereNEDSR*>(NULL));
  ar.register_type(static_cast<Hashed*>(NULL));
  ar.register_type(static_cast<SpaceFilter*>(NULL));
}
