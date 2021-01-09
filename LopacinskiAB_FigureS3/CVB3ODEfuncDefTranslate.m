function [dCdt] = CVB3ODEfuncDefTranslate( t,InitVals,Constants,IFNSwitch,VirResponse,IFNStimulation,IFNStimulationTime )
%Function: Stores the set of ODEs to be run when called from CVB3ODEEval
%
%INPUTS:
    %t: Time
    %InitVals: Vector of simulation initial values
    %Constants: Vector of constants from the ODE solver
    %IFNSwitch: Indicates whether the host-cell anti-viral interferon (IFN) response is enabled. 
    %VirResponse: Indicates whether the viral anti-interferon response is enabled.  
    %IFNStimulation: Indicates if IFN pre-stimulation is simulated.
    %IFNStimulationTime: If IFNStimulation is enabled, indicates the time at which exogenous IFN stimulation occurs. Negative in the case of pre-stimulation.
%    
%OUTPUTS:
    %dCdt: The set of ODEs including time and solutions

%% Imports and assigns all constants and initial values
    
%Assign inputted constants to their respective terms
Constants = num2cell(Constants);
[k_on_DAF,k_off_DAF,k_transloc,k_on_CAR,k_off_CAR,k_internal,InternalVolConv,k_endo_escape, ...
     PolysomeSize,k_T_c_Form,k_Translate,k_P_on,k_P_off,k_N_on,k_N_off,VROVolConv,k_R_Ip_Form,k_N_Transcript,k_P_Transcript,k_R_In_Form, ...
     u_P_cyt,u_P_VRO,u_N_cyt,u_N_VRO,u_T_c,u_VirProt_cyt,u_VirProt_VRO, ...
     k1f_cap,k2f_cap,k1b_cap,k2b_cap,kRNACapBind,kRNACapUnbind,k_Pentamer_on,k_Pentamer_off,u_cap_cyt,u_cap_VRO, ...
     kcat_cleave,Km_cleave,InitRibAvail,u_ISG,ISGformRate,StimISG,ISGBasal,kD_Hill,n_Hill,EC50_Protease,EC50_Translate,OAS_RNAdeg,EC50_RNAdeg,EC50_DetectorDeg,VROFormThreshold] = Constants{:};

%Assign initial conditions to their respective terms
InitVals = num2cell(InitVals);
[uCVB3,uCVB3_Defective,uDAF,bDAF,bDAF_Defective,uDAF_TJ,bDAF_TJ,bDAF_Defective_TJ,uCVB3_TJ,uCVB3_Defective_TJ,uCAR,bCAR,bCAR_Defective,R_p_endo,R_p_endo_Defective, ...
 R_p_cyt,R_p_cyt_Defective,T_c,T_c_Defective,R_n_cyt,R_p_VRO,R_n_VRO,Pol3D_VRO,R_Ip_VRO,R_In_VRO,ATPase2C, ...
 Pentamer_cyt,Pentamer_VRO,RNAPentamer,P2Filled,P3Filled,P4Filled,P5Filled,P6Filled,P7Filled,P8Filled,P9Filled,P10Filled,P11Filled,Virion, ...
 P2Empty,P3Empty,P4Empty,P5Empty,P6Empty,P7Empty,P8Empty,P9Empty,P10Empty,P11Empty,EmptyProvirion, ...
 RibAvail,RibUnavail,Protease,ISGProtein] = InitVals{:};
   
%% Pre-allocates equations for use in the ODE solver

% VRO formation delay
if Pol3D_VRO < VROFormThreshold * 6e-4 %6e-4 is concentration of 1 molecule in nM
    u_P_VRO = u_P_VRO * VROVolConv;
    u_N_VRO = u_N_VRO * VROVolConv;
    u_VirProt_VRO = u_VirProt_VRO * VROVolConv;
    u_cap_VRO = u_cap_VRO * VROVolConv;
    VROVolConv = 1;
end
 
%Receptor binding equations
DAFform = k_on_DAF * uCVB3 * uDAF; %Formation of DAF-CVB3 complex
DAFdiss = k_off_DAF * bDAF; %Dissociation of DAF-CVB3 complex
DAFtoJunc = k_transloc * bDAF; %Movement of DAF-CVB3 complex to the tight junction
DAFoutJunc = k_transloc * uDAF_TJ; %Movement unbound DAF out of the tight junction
DAFform_Junc = k_on_DAF * uCVB3_TJ * uDAF_TJ; %Formation of DAF-CVB3 complex in the tight junction
DAFdiss_Junc = k_off_DAF * bDAF_TJ; %Dissociation of DAF-CVB3 complex in the tight junction
CARform_DAF = k_on_CAR * bDAF_TJ * uCAR; %Transfer of CVB3 from translocalized bDAF onto uCAR
CARform_uCVB3 = k_on_CAR * uCVB3_TJ * uCAR; %Formation of CAR-CVB3 from uCVB3 in the tight junction
CARdiss = k_off_CAR * bCAR; %Dissociation of CAR-CVB3 complex
Internalization = k_internal * bCAR; %Internalization of the virus, leaves uCAR at the membrane
RNARelease = k_endo_escape * R_p_endo; %Release of RNA from internalized virions

%Receptor binding equations involving replication incompetent virions
DAFform_Defect = k_on_DAF * uCVB3_Defective * uDAF; %Formation of DAF-CVB3 complex
DAFdiss_Defect = k_off_DAF * bDAF_Defective; %Dissociation of DAF-CVB3 complex
DAFtoJunc_Defect = k_transloc * bDAF_Defective; %Movement of DAF-CVB3 complex to the tight junction
DAFform_Junc_Defect = k_on_DAF * uCVB3_Defective_TJ * uDAF_TJ; %Formation of DAF-CVB3 complex in the tight junction
DAFdiss_Junc_Defect = k_off_DAF * bDAF_Defective_TJ; %Dissociation of DAF-CVB3 complex in the tight junction
CARform_DAF_Defect = k_on_CAR * bDAF_Defective_TJ * uCAR; %Transfer of CVB3 from translocalized bDAF onto uCAR
CARform_uCVB3_Defect = k_on_CAR * uCVB3_Defective_TJ * uCAR; %Formation of CAR-CVB3 from uCVB3 in the tight junction
CARdiss_Defect = k_off_CAR * bCAR_Defective; %Dissociation of CAR-CVB3 complex
Internalization_Defect = k_internal * bCAR_Defective; %Internalization of the virus, leaves uCAR at the membrane
RNARelease_Defect = k_endo_escape * R_p_endo_Defective; %Release of RNA from internalized virions

%Replication equations
TranslateForm = k_T_c_Form * (RibAvail - T_c - T_c_Defective) * R_p_cyt; %Formation of a translation complex
TranslateDiss = k_Translate * T_c; %Translation completion and complex dissociation
TranslateForm_Defect = k_T_c_Form * (RibAvail - T_c - T_c_Defective) * R_p_cyt_Defective; %Formation of a translation complex
TranslateDiss_Defect = k_Translate * T_c_Defective; %Translation completion and complex dissociation
TranslateDissProt = TranslateDiss; %Separates Protease translation rate for IFN response
PosEnterVRO = k_P_on * R_p_cyt; %%Positive Strand entering the VRO
PosExitVRO = k_P_off * R_p_VRO; %Positive Strand exiting the VRO
MinEnterVRO = k_N_on * R_n_cyt; %Minus Strand entering the VRO
MinExitVRO = k_N_off * R_n_VRO; %Minus Strand exiting the VRO
PosComplexFormVRO = k_R_Ip_Form * R_p_VRO * Pol3D_VRO; %Formation of plus-strand RNA intermediate complexes in the VRO
MinStrandFormVRO = k_N_Transcript * R_Ip_VRO; %Synthesis of negative strands in the VRO
PosStrandFormVRO = k_P_Transcript * R_In_VRO; %Synthesis of positive strands in the VRO
MinComplexFormVRO = k_R_In_Form * R_n_VRO * Pol3D_VRO; %Formation of minus-strand RNA intermediate complexes in the VRO
PosCytDeg = u_P_cyt * R_p_cyt; %Positive Strand degrading in the cytoplasm
PosCytDefectDeg = u_P_cyt * R_p_cyt_Defective; %Defective positive strand degrading in the cytoplasm
PosVRODeg = u_P_VRO * R_p_VRO; %Degradation of positive strands in the VRO
MinCytDeg = u_N_cyt * R_n_cyt; %Minus Strand degrading in the cytoplasm
MinVRODeg = u_N_VRO * R_n_VRO; %Degradation of negative strands in the VRO
TranslateDeg = u_T_c * T_c; %Degradation of translation complexes
TranslateDeg_Defect = u_T_c * T_c_Defective; %Degradation of translation complexes
Pol3DVRODeg = u_VirProt_VRO * Pol3D_VRO; %Degradation of replication enzymes in the VRO
PosComplexDeg = u_P_VRO * R_Ip_VRO; %Degradation of plus-strand RNA intermediate complexes in the VRO
MinComplexDeg = u_N_VRO * R_In_VRO; %Degradation of minus-strand RNA intermediate complexes in the VRO

%Capsid assembly equations
Pent2EmptyAssembly = k1f_cap * Pentamer_VRO^2; %Assembly step for empty capsid molecules (2 unbound pentamers used)
Pent2EmptyDisassembly =  k1b_cap * P2Empty; %Disassembly step for empty capsid molecules in which both unbound pentamers are lost
Pent3EmptyAssembly = k1f_cap * P2Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent3EmptyDisassembly =  k1b_cap * P3Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent4EmptyAssembly = k1f_cap * P3Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent4EmptyDisassembly =  k1b_cap * P4Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent5EmptyAssembly = k1f_cap * P4Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent5EmptyDisassembly =  k1b_cap * P5Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent6EmptyAssembly = k1f_cap * P5Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent6EmptyDisassembly =  k1b_cap * P6Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent7EmptyAssembly = k1f_cap * P6Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent7EmptyDisassembly =  k1b_cap * P7Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent8EmptyAssembly = k1f_cap * P7Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent8EmptyDisassembly =  k1b_cap * P8Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent9EmptyAssembly = k1f_cap * P8Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent9EmptyDisassembly =  k1b_cap * P9Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent10EmptyAssembly = k1f_cap * P9Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent10EmptyDisassembly =  k1b_cap * P10Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
Pent11EmptyAssembly = k1f_cap * P10Empty * Pentamer_VRO; %Assembly step for empty capsid molecules in which an unbound pentamer is added
Pent11EmptyDisassembly =  k1b_cap * P11Empty; %Disassembly step for empty capsid molecules in which an unbound pentamer is lost
EmpProvirionAssembly = k1f_cap * P11Empty * Pentamer_VRO; %Formation of empty provirions from 12 pentamers
EmpProvirionDisassembly = k1b_cap * EmptyProvirion; %Dissociation of empty provirions into 12 pentamers
PentEnterVRO = k_Pentamer_on * Pentamer_cyt * ATPase2C; %Pentamer entering the VRO * VROVolConv
PentExitVRO = k_Pentamer_off * Pentamer_VRO; %Pentamer exiting the VRO
RNAPentBinding = kRNACapBind * Pentamer_VRO * R_p_VRO; %Formation of RNA-Pentamer complexes from species in the VRO
RNAPentUnbinding = kRNACapUnbind * RNAPentamer; %Dissocation of RNA-Pentamer complexes into the cytoplasm
Pent2FilledAssembly = k1f_cap * RNAPentamer * Pentamer_VRO; %Assembly step for filled capsid molecules (1 unbound and 1 RNA-bound pentamer used)
Pent2FilledDisassembly =  k1b_cap * P2Filled; %Disassembly step for filled capsid molecules in which both the unbound and RNA-bound pentamers are lost
Pent3FilledAssembly = k1f_cap * P2Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent3FilledDisassembly =  k1b_cap * P3Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent4FilledAssembly = k1f_cap * P3Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent4FilledDisassembly =  k1b_cap * P4Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent5FilledAssembly = k1f_cap * P4Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent5FilledDisassembly =  k1b_cap * P5Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent6FilledAssembly = k1f_cap * P5Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent6FilledDisassembly =  k1b_cap * P6Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent7FilledAssembly = k1f_cap * P6Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent7FilledDisassembly =  k1b_cap * P7Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent8FilledAssembly = k1f_cap * P7Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent8FilledDisassembly =  k1b_cap * P8Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent9FilledAssembly = k1f_cap * P8Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent9FilledDisassembly =  k1b_cap * P9Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent10FilledAssembly = k1f_cap * P9Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent10FilledDisassembly =  k1b_cap * P10Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
Pent11FilledAssembly = k1f_cap * P10Filled * Pentamer_VRO; %Assembly step for filled capsid molecules in which an unbound pentamer is added 
Pent11FilledDisassembly =  k1b_cap * P11Filled; %Disassembly step for filled capsid molecules in which an unbound pentamer is lost
VirionAssembly = k1f_cap * P11Filled * Pentamer_VRO; %Formation of 12-unit, completed virions from 11-unit filled capsid
VirionDisassembly = 0 * k1b_cap * Virion; %Dissociation of fully encapsidated, RNA-containing virion is assumed to be irreversible
Pent2EmptyFilling = k2f_cap * P2Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent3FilledUnfilling = k2b_cap * P3Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent3EmptyFilling = k2f_cap * P3Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent4FilledUnfilling = k2b_cap * P4Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent4EmptyFilling = k2f_cap * P4Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent5FilledUnfilling = k2b_cap * P5Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent5EmptyFilling = k2f_cap * P5Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent6FilledUnfilling = k2b_cap * P6Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent6EmptyFilling = k2f_cap * P6Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent7FilledUnfilling = k2b_cap * P7Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent7EmptyFilling = k2f_cap * P7Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent8FilledUnfilling = k2b_cap * P8Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent8EmptyFilling = k2f_cap * P8Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent9FilledUnfilling = k2b_cap * P9Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent9EmptyFilling = k2f_cap * P9Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent10FilledUnfilling = k2b_cap * P10Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent10EmptyFilling = k2f_cap * P10Empty * RNAPentamer; %Assembly step for empty capsid in which an RNA-bound pentamer is added
Pent11FilledUnfilling = k2b_cap * P11Filled; %Disassembly step for filled capsid in which an RNA-bound pentamer is lost
Pent11EmptyFilling = k2f_cap * P11Empty * RNAPentamer; %Formation of 12-unit, completed virions from 11-unit empty capsid
PentamerCytDeg = u_cap_cyt * Pentamer_cyt; %Degradation of pentamers in the cytoplasm 
PentamerVRODeg = u_cap_VRO * Pentamer_VRO; %Degradation of pentamers in the VRO 
RNAPentamerDeg = u_cap_VRO * RNAPentamer; %Degradation of RNA-Pentamer complexes 
%All encapsidation intermediates were assumed not to degrade

%EIF4G and PABP cleavage equations
RiboChange = kcat_cleave * Protease * RibUnavail * PolysomeSize / (RibUnavail * PolysomeSize + Km_cleave); %Michaelis-Menten equation for protease activity against eIF4G 
ProtCytDeg = u_VirProt_cyt * Protease; %Degradation of viral proteases in the cytoplasm

%IFN-B response equations
TotalDoubleStrands = (R_In_VRO + R_Ip_VRO); %Summed for subsequent use in Viral detection equation

%Set VirDetection (dsRNA sensing mediated IFN response) and StimISG (exogenous IFN stimulation) depending on the amount of dsRNA and the simulation time respectively.
if strcmpi(IFNSwitch,'on') %If the interferon response module is active
    
    VirDetection = ((TotalDoubleStrands)^n_Hill)/(kD_Hill^n_Hill + (TotalDoubleStrands)^n_Hill); %VirDetection is defined by a Hill equation for detection of dsRNA
    
    if strcmpi(IFNStimulation,'on') %If there is exogenous interferon stimulation
        if IFNStimulationTime <= 0 % and if the stimulation time occurs prior to infection (pre-stimulation)
            StimISG = 1; % then StimISG is always maximally on.
            VirDetection = 0; %dsRNA mediated IFN response is superseded by the exogenous response.
        elseif IFNStimulationTime > 0 && t < IFNStimulationTime %If the time of stimulation has not occured yet,
            StimISG = 0; %then StimISG is off.
        elseif IFNStimulationTime > 0 && t >= IFNStimulationTime %If the time of stimulation has occured, 
            StimISG = 1; %then StimISG is maximally on.
            VirDetection = 0;%dsRNA mediated IFN response is superseded by the exogenous response.
        end
    end
    
    if strcmpi(VirResponse,'on')
        VirDetection = VirDetection*(1 - Protease/(Protease + EC50_DetectorDeg)); %Scales down viral detection for IFN response if the virus's response to the IFN response is active.
    end    
else
    VirDetection = 0; %Sets detection to 0 in the event of no IFN response
end

ISGform = (VirDetection + StimISG) * ISGformRate * RibUnavail; %Formation of ISGs can be induced by either VirDetection or exogenous StimISG.
ISGDegradation = u_ISG * ISGProtein; %Degradation of ISG proteins in the cytoplasm
ISGBasalRate = u_ISG * ISGBasal; %Maintains the basal level of ISG proteins against degradation
if strcmpi(IFNSwitch,'on') %Impedes virus processes if IFN response is turned on
    TranslateForm = TranslateForm * (1 - ISGProtein/(EC50_Translate + ISGProtein)); %Formation of translation complexes inhibited by interferon
    TranslateDissProt = TranslateDissProt * (1 - ISGProtein/(EC50_Protease + ISGProtein)); %ISG15 deactivation of proteases by ISGylation
    
    PosCytDeg = PosCytDeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
    PosCytDefectDeg = PosCytDeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
    PosVRODeg = PosVRODeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
    MinCytDeg = MinCytDeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
    MinVRODeg = MinVRODeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
    TranslateDeg = TranslateDeg * (1 + (OAS_RNAdeg-1)*ISGProtein/(EC50_RNAdeg + ISGProtein));   %Ramp degradation up 5-fold higher as ISGs are expressed, modeling OASs
end


%% Defining the system of ODEs

dCdt = [DAFdiss - DAFform; %uCVB3
    DAFdiss_Defect - DAFform_Defect; %uCVB3_Defective
    DAFdiss - DAFform - DAFform_Defect + DAFdiss_Defect + DAFoutJunc; %uDAF
    DAFform - DAFdiss - DAFtoJunc; %bDAF
    DAFform_Defect - DAFdiss_Defect - DAFtoJunc_Defect; %bDAF_Defective
    CARform_DAF - DAFoutJunc + DAFdiss_Junc - DAFform_Junc + CARform_DAF_Defect + DAFdiss_Junc_Defect - DAFform_Junc_Defect; %uDAF_TJ
    DAFtoJunc + DAFform_Junc - DAFdiss_Junc - CARform_DAF; %bDAF_TJ
    DAFtoJunc_Defect + DAFform_Junc_Defect - DAFdiss_Junc_Defect - CARform_DAF_Defect; %bDAF_Defective_TJ
    DAFdiss_Junc + CARdiss - CARform_uCVB3 - DAFform_Junc; %uCVB3_TJ
    DAFdiss_Junc_Defect + CARdiss_Defect - CARform_uCVB3_Defect - DAFform_Junc_Defect; %uCVB3_Defective_TJ
    CARdiss - CARform_DAF - CARform_uCVB3 + CARdiss_Defect - CARform_DAF_Defect - CARform_uCVB3_Defect; %uCAR
    CARform_DAF + CARform_uCVB3 - CARdiss - Internalization; %bCAR
    CARform_DAF_Defect + CARform_uCVB3_Defect - CARdiss_Defect - Internalization_Defect; %bCAR_Defective
    InternalVolConv * Internalization - RNARelease; %R_p_endo
    InternalVolConv * Internalization_Defect - RNARelease_Defect; %R_p_endo_Defective
    RNARelease + PosExitVRO / VROVolConv + RNAPentUnbinding / VROVolConv - TranslateForm - PosEnterVRO - PosCytDeg; %R_p_cyt
    RNARelease_Defect - TranslateForm_Defect - PosCytDefectDeg + TranslateDiss_Defect/PolysomeSize; %R_p_cyt_Defective
    TranslateForm - TranslateDiss/PolysomeSize - TranslateDeg; %T_c
    TranslateForm_Defect - TranslateDiss_Defect/PolysomeSize - TranslateDeg_Defect; %T_c_Defect
    MinExitVRO / VROVolConv - MinEnterVRO - MinCytDeg; %R_n_cyt
    TranslateDiss/PolysomeSize * VROVolConv + PosStrandFormVRO + MinStrandFormVRO - PosComplexFormVRO - RNAPentBinding + PosEnterVRO * VROVolConv - PosExitVRO - PosVRODeg; %R_p_VRO
    MinStrandFormVRO + MinEnterVRO * VROVolConv - MinExitVRO - MinComplexFormVRO - MinVRODeg; %R_n_VRO
    TranslateDiss * VROVolConv + MinStrandFormVRO - PosComplexFormVRO - MinComplexFormVRO + PosComplexDeg + MinComplexDeg - Pol3DVRODeg; %Pol3D_VRO
    PosComplexFormVRO - MinStrandFormVRO - PosComplexDeg; %R_Ip_VRO
    MinComplexFormVRO - MinComplexDeg; %R_In_VRO
    4/5*TranslateDiss * VROVolConv + PentExitVRO - PentEnterVRO * VROVolConv + RNAPentBinding ...
        + Pent2FilledAssembly + Pent2EmptyAssembly + Pent3FilledAssembly + Pent3EmptyAssembly + Pent4FilledAssembly + Pent4EmptyAssembly + Pent5FilledAssembly + Pent5EmptyAssembly ...
        + Pent6FilledAssembly + Pent6EmptyAssembly + Pent7FilledAssembly + Pent7EmptyAssembly + Pent8FilledAssembly + Pent8EmptyAssembly + Pent9FilledAssembly + Pent9EmptyAssembly ...
        + Pent10FilledAssembly + Pent10EmptyAssembly + Pent11FilledAssembly + Pent11EmptyAssembly + VirionAssembly + EmpProvirionAssembly ...
        + Pent2EmptyFilling + Pent3EmptyFilling + Pent4EmptyFilling + Pent5EmptyFilling + Pent6EmptyFilling + Pent7EmptyFilling + Pent8EmptyFilling + Pent9EmptyFilling + Pent10EmptyFilling + Pent11EmptyFilling ...
        - 1/3*Pent3FilledUnfilling - 1/4*Pent4FilledUnfilling - 1/5*Pent5FilledUnfilling - 1/6*Pent6FilledUnfilling - 1/7*Pent7FilledUnfilling - 1/8*Pent8FilledUnfilling - 1/9*Pent9FilledUnfilling ...
        - 1/10*Pent10FilledUnfilling - 1/11*Pent11FilledUnfilling; %ATPase2C (4/5 because 1 is made per polyprotein, but 1/5 is allocated as a Pentamer_VRO
    RNAPentUnbinding / VROVolConv + PentExitVRO / VROVolConv - PentEnterVRO - PentamerCytDeg ...
        + (Pent2FilledDisassembly + Pent2EmptyDisassembly + 2/3*Pent3FilledDisassembly + Pent3EmptyDisassembly + 3/4*Pent4FilledDisassembly + Pent4EmptyDisassembly + 4/5*Pent5FilledDisassembly + Pent5EmptyDisassembly ...
        + 5/6*Pent6FilledDisassembly + Pent6EmptyDisassembly + 6/7*Pent7FilledDisassembly + Pent7EmptyDisassembly + 7/8*Pent8FilledDisassembly + Pent8EmptyDisassembly + 8/9*Pent9FilledDisassembly + Pent9EmptyDisassembly ...
        + 9/10*Pent10FilledDisassembly + Pent10EmptyDisassembly + 10/11*Pent11FilledDisassembly + Pent11EmptyDisassembly + VirionDisassembly + EmpProvirionDisassembly)/VROVolConv; %Pentamer_cyt
    TranslateDiss/5 * VROVolConv + PentEnterVRO * VROVolConv - PentExitVRO - RNAPentBinding - PentamerVRODeg + Pent2EmptyDisassembly ...
        - Pent2FilledAssembly - 2*Pent2EmptyAssembly - Pent3FilledAssembly - Pent3EmptyAssembly - Pent4FilledAssembly - Pent4EmptyAssembly - Pent5FilledAssembly - Pent5EmptyAssembly ...
        - Pent6FilledAssembly - Pent6EmptyAssembly - Pent7FilledAssembly - Pent7EmptyAssembly - Pent8FilledAssembly - Pent8EmptyAssembly - Pent9FilledAssembly - Pent9EmptyAssembly ...
        - Pent10FilledAssembly - Pent10EmptyAssembly - Pent11FilledAssembly - Pent11EmptyAssembly - VirionAssembly - EmpProvirionAssembly; %Pentamer_VRO = Pentamer:ATPase2C
    RNAPentBinding - RNAPentUnbinding - RNAPentamerDeg - Pent2FilledAssembly + Pent2FilledDisassembly - Pent2EmptyFilling ...
        + 1/3*Pent3FilledUnfilling - Pent3EmptyFilling + 1/4*Pent4FilledUnfilling - Pent4EmptyFilling + 1/5*Pent5FilledUnfilling - Pent5EmptyFilling + 1/6*Pent6FilledUnfilling - Pent6EmptyFilling ...
        + 1/7*Pent7FilledUnfilling - Pent7EmptyFilling + 1/8*Pent8FilledUnfilling - Pent8EmptyFilling + 1/9*Pent9FilledUnfilling - Pent9EmptyFilling + 1/10*Pent10FilledUnfilling - Pent10EmptyFilling ...
        + 1/11*Pent11FilledUnfilling - Pent11EmptyFilling; %RNAPentamer = Pentamer:R_p_VRO
    Pent2FilledAssembly - Pent2FilledDisassembly + 2/3*Pent3FilledDisassembly - Pent3FilledAssembly; %P2Filled = 2xPentamer:R_p_VRO
    Pent3FilledAssembly - 2/3*Pent3FilledDisassembly + 3/4*Pent4FilledDisassembly - Pent4FilledAssembly - 1/3*Pent3FilledUnfilling + Pent2EmptyFilling; %P3Filled = 3xPentamer:R_p_VRO
    Pent4FilledAssembly - 3/4*Pent4FilledDisassembly + 4/5*Pent5FilledDisassembly - Pent5FilledAssembly - 1/4*Pent4FilledUnfilling + Pent3EmptyFilling; %P4Filled = 4xPentamer:R_p_VRO
    Pent5FilledAssembly - 4/5*Pent5FilledDisassembly + 5/6*Pent6FilledDisassembly - Pent6FilledAssembly - 1/5*Pent5FilledUnfilling + Pent4EmptyFilling; %P5Filled = 5xPentamer:R_p_VRO
    Pent6FilledAssembly - 5/6*Pent6FilledDisassembly + 6/7*Pent7FilledDisassembly - Pent7FilledAssembly - 1/6*Pent6FilledUnfilling + Pent5EmptyFilling; %P6Filled = 6xPentamer:R_p_VRO
    Pent7FilledAssembly - 6/7*Pent7FilledDisassembly + 7/8*Pent8FilledDisassembly - Pent8FilledAssembly - 1/7*Pent7FilledUnfilling + Pent6EmptyFilling; %P7Filled = 7xPentamer:R_p_VRO
    Pent8FilledAssembly - 7/8*Pent8FilledDisassembly + 8/9*Pent9FilledDisassembly - Pent9FilledAssembly - 1/8*Pent8FilledUnfilling + Pent7EmptyFilling; %P8Filled = 8xPentamer:R_p_VRO
    Pent9FilledAssembly - 8/9*Pent9FilledDisassembly + 9/10*Pent10FilledDisassembly - Pent10FilledAssembly - 1/9*Pent9FilledUnfilling + Pent8EmptyFilling; %P9Filled = 9xPentamer:R_p_VRO
    Pent10FilledAssembly - 9/10*Pent10FilledDisassembly + 10/11*Pent11FilledDisassembly - Pent11FilledAssembly - 1/10*Pent10FilledUnfilling + Pent9EmptyFilling; %P10Filled = 10xPentamer:R_p_VRO
    Pent11FilledAssembly - 10/11*Pent11FilledDisassembly + VirionDisassembly - VirionAssembly - 1/11*Pent11FilledUnfilling + Pent10EmptyFilling; %P11Filled = 11xPentamer:R_p_VRO
    VirionAssembly - VirionDisassembly + Pent11EmptyFilling; %Virion = 12xPentamer:R_p_VRO
    Pent2EmptyAssembly - Pent2EmptyDisassembly + Pent3EmptyDisassembly - Pent3EmptyAssembly + 1/3*Pent3FilledUnfilling - Pent2EmptyFilling; %P2Empty = 2xPentamer:ATPase2C
    Pent3EmptyAssembly - Pent3EmptyDisassembly + Pent4EmptyDisassembly - Pent4EmptyAssembly + 1/4*Pent4FilledUnfilling - Pent3EmptyFilling; %P3Empty = 3xPentamer:ATPase2C
    Pent4EmptyAssembly - Pent4EmptyDisassembly + Pent5EmptyDisassembly - Pent5EmptyAssembly + 1/5*Pent5FilledUnfilling - Pent4EmptyFilling; %P4Empty = 4xPentamer:ATPase2C
    Pent5EmptyAssembly - Pent5EmptyDisassembly + Pent6EmptyDisassembly - Pent6EmptyAssembly + 1/6*Pent6FilledUnfilling - Pent5EmptyFilling; %P5Empty = 5xPentamer:ATPase2C
    Pent6EmptyAssembly - Pent6EmptyDisassembly + Pent7EmptyDisassembly - Pent7EmptyAssembly + 1/7*Pent7FilledUnfilling - Pent6EmptyFilling; %P6Empty = 6xPentamer:ATPase2C
    Pent7EmptyAssembly - Pent7EmptyDisassembly + Pent8EmptyDisassembly - Pent8EmptyAssembly + 1/8*Pent8FilledUnfilling - Pent7EmptyFilling; %P7Empty = 7xPentamer:ATPase2C
    Pent8EmptyAssembly - Pent8EmptyDisassembly + Pent9EmptyDisassembly - Pent9EmptyAssembly + 1/9*Pent9FilledUnfilling - Pent8EmptyFilling; %P8Empty = 8xPentamer:ATPase2C
    Pent9EmptyAssembly - Pent9EmptyDisassembly + Pent10EmptyDisassembly - Pent10EmptyAssembly + 1/10*Pent10FilledUnfilling - Pent9EmptyFilling; %P9Empty = 9xPentamer:ATPase2C
    Pent10EmptyAssembly - Pent10EmptyDisassembly + Pent11EmptyDisassembly - Pent11EmptyAssembly + 1/11*Pent11FilledUnfilling - Pent10EmptyFilling; %P10Empty = 10xPentamer:ATPase2C
    Pent11EmptyAssembly - Pent11EmptyDisassembly + EmpProvirionDisassembly - EmpProvirionAssembly - Pent11EmptyFilling; %P11Empty = 11xPentamer:ATPase2C   
    EmpProvirionAssembly - EmpProvirionDisassembly; %EmptyProvirion = 12xPentamer:ATPase2C
    RiboChange; %RibAvail
    -RiboChange; %RibUnavail
    TranslateDissProt - ProtCytDeg; %Protease
    ISGform + ISGBasalRate - ISGDegradation]; %ISGProtein

end