# Device and service model configuration for the UC, MD, and ED PSI templates.
# Consumed by create_uc_template, create_md_template, and create_ed_template in
# markets_simulation/siip_simulation_definition.jl.

const UC_DEVICE_MODELS = Dict(
    PSY.ThermalStandard => PSI.ThermalBasicUnitCommitment,
    ThermalFastStartSIIP => PSI.ThermalBasicUnitCommitment,
    PSY.ThermalMultiStart => PSI.ThermalBasicUnitCommitment,
    PSY.RenewableDispatch => PSI.RenewableFullDispatch,
    PSY.RenewableNonDispatch => PSI.FixedOutput,
    PSY.StandardLoad => PSI.StaticPowerLoad,
    PSY.HydroTurbine => HSI.HydroCommitmentRunOfRiver,
    PSY.HydroDispatch => HSI.HydroCommitmentRunOfRiver, # TODO: check which hydro device we have
    PSY.EnergyReservoirStorage => SSI.StorageDispatchWithReserves,
    PSY.Line => PSI.StaticBranch,
    PSY.Transformer2W => PSI.StaticBranch,
    PSY.TapTransformer => PSI.StaticBranch,
    PSY.TwoTerminalGenericHVDCLine => PSI.HVDCTwoTerminalLossless,
)

# NOTE: unlike UC_DEVICE_MODELS, this has no PSY.ThermalMultiStart entry — matches
# the pre-existing behavior of create_md_template, not necessarily intentional.
const MD_DEVICE_MODELS = Dict(
    PSY.ThermalStandard => PSI.ThermalBasicUnitCommitment,
    ThermalFastStartSIIP => PSI.ThermalBasicUnitCommitment,
    PSY.RenewableDispatch => PSI.RenewableFullDispatch,
    PSY.RenewableNonDispatch => PSI.FixedOutput,
    PSY.StandardLoad => PSI.StaticPowerLoad,
    PSY.HydroTurbine => HSI.HydroCommitmentRunOfRiver,
    PSY.HydroDispatch => HSI.HydroCommitmentRunOfRiver, # TODO: check which hydro device we have
    PSY.EnergyReservoirStorage => SSI.StorageDispatchWithReserves,
    PSY.Line => PSI.StaticBranch,
    PSY.Transformer2W => PSI.StaticBranch,
    PSY.TapTransformer => PSI.StaticBranch,
    PSY.TwoTerminalGenericHVDCLine => PSI.HVDCTwoTerminalLossless,
)

const ED_DEVICE_MODELS = Dict(
    PSY.ThermalStandard => PSI.ThermalBasicDispatch,
    ThermalFastStartSIIP => PSI.ThermalBasicUnitCommitment,
    PSY.ThermalMultiStart => PSI.ThermalBasicUnitCommitment,
    PSY.RenewableDispatch => PSI.RenewableFullDispatch,
    PSY.RenewableNonDispatch => PSI.FixedOutput,
    PSY.StandardLoad => PSI.StaticPowerLoad,
    PSY.HydroTurbine => HSI.HydroDispatchRunOfRiver,
    PSY.HydroDispatch => HSI.HydroDispatchRunOfRiver, # TODO: check which hydro device we have
    PSY.EnergyReservoirStorage => SSI.StorageDispatchWithReserves,
    PSY.Line => PSI.StaticBranch,
    PSY.Transformer2W => PSI.StaticBranch,
    PSY.TapTransformer => PSI.StaticBranch,
    PSY.TwoTerminalGenericHVDCLine => PSI.HVDCTwoTerminalLossless,
)

const REG_UP_SERVICE_MODEL = (
    component = PSY.VariableReserve{PSY.ReserveUp},
    formulation = PSI.RangeReserve,
    name = "Reg_Up",
    use_slacks = true,
    duals = [PSI.RequirementConstraint],
)

const REG_DOWN_SERVICE_MODEL = (
    component = PSY.VariableReserve{PSY.ReserveDown},
    formulation = PSI.RangeReserve,
    name = "Reg_Down",
    use_slacks = true,
    duals = [PSI.RequirementConstraint],
)

const ORDC_SYNCHRONOUS_SERVICE_MODEL = (
    component = PSY.ReserveDemandCurve{PSY.ReserveUp},
    formulation = PSI.StepwiseCostReserve,
    name = "Synchronous",
    use_slacks = true,
    duals = [PSI.RequirementConstraint],
)

const ORDC_PRIMARY_SERVICE_MODEL = (
    component = PSY.ReserveDemandCurve{PSY.ReserveUp},
    formulation = PSI.StepwiseCostReserve,
    name = "Primary",
    use_slacks = true,
    duals = [PSI.RequirementConstraint],
)

const UC_SERVICE_MODELS = [
    REG_UP_SERVICE_MODEL,
    REG_DOWN_SERVICE_MODEL,
    ORDC_SYNCHRONOUS_SERVICE_MODEL,
    ORDC_PRIMARY_SERVICE_MODEL,
]

const MD_SERVICE_MODELS = UC_SERVICE_MODELS

# ED always gets Reg_Up/Reg_Down. The ORDC service models are added only when
# inertia_product is non-empty (see create_ed_template).
const ED_SERVICE_MODELS_BASE = [REG_UP_SERVICE_MODEL, REG_DOWN_SERVICE_MODEL]
const ED_SERVICE_MODELS_ORDC = [ORDC_SYNCHRONOUS_SERVICE_MODEL, ORDC_PRIMARY_SERVICE_MODEL]

function apply_device_models!(template::PSI.ProblemTemplate, device_models::Dict)
    for (device_type, formulation) in device_models
        PSI.set_device_model!(template, device_type, formulation)
    end
    return
end

function apply_service_models!(template::PSI.ProblemTemplate, service_models::Vector)
    for svc in service_models
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                svc.component,
                svc.formulation,
                svc.name;
                use_slacks = svc.use_slacks,
                duals = svc.duals,
            ),
        )
    end
    return
end
