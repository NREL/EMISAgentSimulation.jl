function convert(::Type{Project{Existing}}, gen::ThermalGenEMIS{BuildPhase})
    return ThermalGenEMIS{Existing}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::ThermalGenEMIS{BuildPhase})
    return ThermalGenEMIS{Planned}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Queue}}, gen::ThermalGenEMIS{BuildPhase})
    return ThermalGenEMIS{Queue}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Option}}, gen::ThermalGenEMIS{BuildPhase})
    return ThermalGenEMIS{Option}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::ThermalGenEMIS{BuildPhase})
    return ThermalGenEMIS{Retired}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Existing}}, gen::RenewableGenEMIS{BuildPhase})
    return RenewableGenEMIS{Existing}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::RenewableGenEMIS{BuildPhase})
    return RenewableGenEMIS{Planned}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Queue}}, gen::RenewableGenEMIS{BuildPhase})
    return RenewableGenEMIS{Queue}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Option}}, gen::RenewableGenEMIS{BuildPhase})
    return RenewableGenEMIS{Option}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::RenewableGenEMIS{BuildPhase})
    return RenewableGenEMIS{Retired}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Existing}}, gen::HydroGenEMIS{BuildPhase})
    return HydroGenEMIS{Existing}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::HydroGenEMIS{BuildPhase})
    return HydroGenEMIS{Planned}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Queue}}, gen::HydroGenEMIS{BuildPhase})
    return HydroGenEMIS{Queue}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Option}}, gen::HydroGenEMIS{BuildPhase})
    return HydroGenEMIS{Option}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::HydroGenEMIS{BuildPhase})
    return HydroGenEMIS{Retired}(gen.name,
                               gen.tech,
                               gen.decision_year,
                               gen.construction_year,
                               gen.retirement_year,
                               gen.end_life_year,
                               gen.products,
                               gen.finance_data)
end

function convert(::Type{Project{Existing}}, battery::BatteryEMIS{BuildPhase})
    return BatteryEMIS{Existing}(battery.name,
                               battery.tech,
                               battery.decision_year,
                               battery.construction_year,
                               battery.retirement_year,
                               battery.end_life_year,
                               battery.products,
                               battery.finance_data)
end

function convert(::Type{Project{Planned}}, battery::BatteryEMIS{BuildPhase})
    return BatteryEMIS{Planned}(battery.name,
                               battery.tech,
                               battery.decision_year,
                               battery.construction_year,
                               battery.retirement_year,
                               battery.end_life_year,
                               battery.products,
                               battery.finance_data)
end

function convert(::Type{Project{Queue}}, battery::BatteryEMIS{BuildPhase})
    return BatteryEMIS{Queue}(battery.name,
                               battery.tech,
                               battery.decision_year,
                               battery.construction_year,
                               battery.retirement_year,
                               battery.end_life_year,
                               battery.products,
                               battery.finance_data)
end

function convert(::Type{Project{Option}}, battery::BatteryEMIS{BuildPhase})
    return BatteryEMIS{Option}(battery.name,
                               battery.tech,
                               battery.decision_year,
                               battery.construction_year,
                               battery.retirement_year,
                               battery.end_life_year,
                               battery.products,
                               battery.finance_data)
end

function convert(::Type{Project{Retired}}, battery::BatteryEMIS{BuildPhase})
    return BatteryEMIS{Retired}(battery.name,
                               battery.tech,
                               battery.decision_year,
                               battery.construction_year,
                               battery.retirement_year,
                               battery.end_life_year,
                               battery.products,
                               battery.finance_data)
end


function convert(::Type{Project{Queue}}, gen::ThermalGenEMIS{Option})
    return ThermalGenEMIS{Queue}(gen.name,
                            gen.tech,
                            gen.decision_year,
                            gen.construction_year,
                            gen.retirement_year,
                            gen.end_life_year,
                            gen.products,
                            gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::ThermalGenEMIS{Queue})
    return ThermalGenEMIS{Planned}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Existing}}, gen::ThermalGenEMIS{Planned})
    return ThermalGenEMIS{Existing}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end


function convert(::Type{Project{Retired}}, gen::ThermalGenEMIS{Existing})
    return ThermalGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::ThermalGenEMIS{Planned})
    return ThermalGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::ThermalGenEMIS{Queue})
    return ThermalGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Queue}}, gen::RenewableGenEMIS{Option})
    return RenewableGenEMIS{Queue}(gen.name,
                            gen.tech,
                            gen.decision_year,
                            gen.construction_year,
                            gen.retirement_year,
                            gen.end_life_year,
                            gen.products,
                            gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::RenewableGenEMIS{Queue})
    return RenewableGenEMIS{Planned}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Existing}}, gen::RenewableGenEMIS{Planned})
    return RenewableGenEMIS{Existing}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end


function convert(::Type{Project{Retired}}, gen::RenewableGenEMIS{Existing})
    return RenewableGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::RenewableGenEMIS{Planned})
    return RenewableGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::RenewableGenEMIS{Queue})
    return RenewableGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Queue}}, gen::HydroGenEMIS{Option})
    return HydroGenEMIS{Queue}(gen.name,
                            gen.tech,
                            gen.decision_year,
                            gen.construction_year,
                            gen.retirement_year,
                            gen.end_life_year,
                            gen.products,
                            gen.finance_data)
end

function convert(::Type{Project{Planned}}, gen::HydroGenEMIS{Queue})
    return HydroGenEMIS{Planned}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Existing}}, gen::HydroGenEMIS{Planned})
    return HydroGenEMIS{Existing}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::HydroGenEMIS{Existing})
    return HydroGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::HydroGenEMIS{Planned})
    return HydroGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Retired}}, gen::HydroGenEMIS{Queue})
    return HydroGenEMIS{Retired}(gen.name,
                              gen.tech,
                              gen.decision_year,
                              gen.construction_year,
                              gen.retirement_year,
                              gen.end_life_year,
                              gen.products,
                              gen.finance_data)
end

function convert(::Type{Project{Queue}}, battery::BatteryEMIS{Option})
    return BatteryEMIS{Queue}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end

function convert(::Type{Project{Planned}}, battery::BatteryEMIS{Queue})
    return BatteryEMIS{Planned}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end

function convert(::Type{Project{Existing}}, battery::BatteryEMIS{Planned})
    return BatteryEMIS{Existing}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end


function convert(::Type{Project{Retired}}, battery::BatteryEMIS{Existing})
    return BatteryEMIS{Retired}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end

function convert(::Type{Project{Retired}}, battery::BatteryEMIS{Planned})
    return BatteryEMIS{Retired}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end

function convert(::Type{Project{Retired}}, battery::BatteryEMIS{Queue})
    return BatteryEMIS{Retired}(battery.name,
                            battery.tech,
                            battery.decision_year,
                            battery.construction_year,
                            battery.retirement_year,
                            battery.end_life_year,
                            battery.products,
                            battery.finance_data)
end
