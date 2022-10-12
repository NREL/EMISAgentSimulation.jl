abstract type BuildPhase end

abstract type Existing <: BuildPhase end    #Existing projects.
abstract type Planned <: BuildPhase end     #Projects under construction.
abstract type Queue <: BuildPhase end       #Projects in queue.
abstract type Option <: BuildPhase end      #Projects which can be built.
abstract type Retired <: BuildPhase end     #Retired projects.
