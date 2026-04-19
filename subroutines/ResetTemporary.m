function ResetTemporary()
    
    global time_pcg time_Psis time_Psi time_Phis time_Phi

    time_Psis = [];
    time_Psi  = [];
    time_Phis = [];
    time_Phi  = [];
    time_pcg.invtildeM = EventTime.pcg();
    time_pcg.invtildeN = EventTime.pcg();
end