package edu.scripps.yates.nucleome.turboID;

public enum IP_Channel_Ori {
	ori_T1(IPExperimentType.NU), ori_T2(null), ori_T3(IPExperimentType.BKGD300), ori_T4(IPExperimentType.EMD150a),
	ori_T5(IPExperimentType.EMD300a), ori_T6(IPExperimentType.EMD150b), ori_T7(IPExperimentType.EMD300b),
	ori_T8(IPExperimentType.Bt300)
//	,	ori_T9(IPExperimentType.NU), ori_T10(IPExperimentType.NU), ori_T11(IPExperimentType.NU)
	;

	private final IPExperimentType expType;

	private IP_Channel_Ori(IPExperimentType expType) {
		this.expType = expType;
	}

	public IPExperimentType getExpType() {
		return expType;
	}

}
