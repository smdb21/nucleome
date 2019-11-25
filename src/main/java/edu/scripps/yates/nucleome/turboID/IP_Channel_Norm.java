package edu.scripps.yates.nucleome.turboID;

public enum IP_Channel_Norm {
	norm_T1(IPExperimentType.NU), norm_T2(null), norm_T3(IPExperimentType.BKGD300), norm_T4(IPExperimentType.EMD150a),
	norm_T5(IPExperimentType.EMD300a), norm_T6(IPExperimentType.EMD150b), norm_T7(IPExperimentType.EMD300b),
	norm_T8(IPExperimentType.Bt300),
//			norm_T9(IPExperimentType., norm_T10(IPExperimentType., norm_T11(IPExperimentType.;
	;
	private final IPExperimentType expType;

	private IP_Channel_Norm(IPExperimentType expType) {
		this.expType = expType;
	}

	public IPExperimentType getExpType() {
		return expType;
	}
}
